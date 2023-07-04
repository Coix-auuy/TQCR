import gurobipy as gp
from common.data_reader import DataReader


class MIQPSolver:
    def solve(self, file_path):
        """
        求解
        :param file_path: 求解路径
        :return:
        """
        data = DataReader(file_path)
        t_num = data.t_num
        u_num = data.u_num
        load = data.load
        reserve = data.reserve
        r_up = data.r_up
        r_down = data.r_down
        a = data.a
        b = data.b
        c = data.c
        min_output = data.min_output
        max_output = data.max_output
        t_init = data.t_init
        p_init = data.p_init
        t_on = data.t_on
        t_off = data.t_off
        c_h = data.c_h
        c_c = data.c_c
        t_c = data.t_c
        r_on = data.r_on
        r_off = data.r_off
        u0 = data.u0

        try:
            # Model start
            m = gp.Model("MIQP")
            # Create variables
            rc = list([(i, t) for i in range(1, u_num + 1) for t in range(1, t_num + 1)])
            u = m.addVars(rc, vtype=gp.GRB.BINARY, name="u")
            v = m.addVars(rc, vtype=gp.GRB.BINARY, name="v")
            w = m.addVars(rc, vtype=gp.GRB.BINARY, name="w")
            p = m.addVars(rc, vtype=gp.GRB.CONTINUOUS, name="p")
            sc = m.addVars(rc, vtype=gp.GRB.CONTINUOUS, name="sc")

            # Objective function (1)
            m.setObjective(gp.quicksum(
                a[i - 1] * u[i, t] + b[i - 1] * p[i, t] + c[i - 1] * (p[i, t] ** 2) + sc[i, t] for i in range(1, u_num + 1) for t in
                range(1, t_num + 1)), gp.GRB.MINIMIZE)

            # Constraints
            # startup cost (2)(3)
            m.addConstrs((sc[i, t] >= v[i, t] * c_h[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                         name="startup cost constraints 1")
            m.addConstrs(
                (sc[i, t] >= c_c[i - 1] * (v[i, t] - gp.quicksum(w[i, k] for k in range(max(t - t_off[i - 1] - t_c[i - 1], 1), t)) - (
                    1 if t - t_off[i - 1] - t_c[i - 1] <= 0 and max(0, -t_init[i - 1]) < abs(
                        t - t_off[i - 1] - t_c[i - 1] - 1) + 1 else 0))
                 for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="startup cost constraints 2")

            # power balance constraints (4)
            m.addConstrs((gp.quicksum(p[i, t] for i in range(1, u_num + 1)) == load[t - 1] for t in range(1, t_num + 1)),
                         name="power balance constraints")

            # system spinning reserve requirements (5)
            m.addConstrs((
                gp.quicksum(u[i, t] * max_output[i - 1] for i in range(1, u_num + 1)) >= load[t - 1] +
                reserve[t - 1] for t in range(1, t_num + 1)), name="spinning reserve requirements")

            # unit generation limits (6)
            m.addConstrs((p[i, t] >= u[i, t] * min_output[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                         name="unit generation limits 1")

            # unit generation limits (6)
            m.addConstrs((p[i, t] <= u[i, t] * max_output[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                         name="unit generation limits 2")

            # ramp rate limits (7) (8)
            m.addConstrs(
                (p[i, t] - (p[i, t - 1] if t > 1 else p_init[i - 1]) <= (u[i, t - 1] if t > 1 else u0[i - 1]) * r_up[i - 1] + v[i, t] *
                 r_on[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="ramp rate limits 1")
            m.addConstrs(
                ((p[i, t - 1] if t > 1 else p_init[i - 1]) - p[i, t] <= u[i, t] * r_down[i - 1] + w[i, t] * r_off[i - 1] for
                 i in
                 range(1, u_num + 1) for t in range(1, t_num + 1)), name="ramp rate limits 2")

            # minimum up and downtime constraints can be formulated (9)(10)
            W = []
            for i in range(1, u_num + 1):
                W.append(max(
                    0, min(t_num, u0[i - 1] * (t_on[i - 1] - t_init[i - 1]))))
            L = []
            for i in range(1, u_num + 1):
                L.append(max(0, min(
                    t_num, (1 - u0[i - 1]) * (t_off[i - 1] + t_init[i - 1]))))
            m.addConstrs(
                (gp.quicksum(v[i, j] for j in range(max(0, t - t_on[i - 1]) + 1, t + 1)) <= u[i, t] for i in range(1, u_num + 1) for t
                 in range(W[i - 1] + 1, t_num + 1)), name=" minimum on time constraints")
            m.addConstrs(
                (gp.quicksum(w[i, j] for j in range(max(0, t - t_off[i - 1]) + 1, t + 1)) <= 1 - u[i, t] for i in range(1, u_num + 1)
                 for t in range(L[i - 1] + 1, t_num + 1)), name="maximum off time constraints")

            #  state variables and logical constraints (11)
            m.addConstrs(
                (v[i, t] - w[i, t] == u[i, t] - (u[i, t - 1] if t > 1 else u0[i - 1]) for i in range(1, u_num + 1) for t in
                 range(1, t_num + 1)),
                name="logical constraints")

            # initial status of units
            m.addConstrs(
                (u[i, t] == u0[i - 1] for i in range(1, u_num + 1)
                 for t in range(1, W[i - 1] + L[i - 1] + 1)),
                name="initial status constraints")
            # Model end

            # Compute optimal cost
            # 设置 gurobi 参数
            for key, value in self.gurobi_params.items():
                m.setParam(key, value)
            m.update()
            m.optimize()
            # 获取求解后的结果
            self.ObjVal = m.getAttr("ObjVal")
            self.MIPGap = m.getAttr("MIPGap")
            self.Runtime = m.getAttr("Runtime")
            self.p = [[0] * (t_num + 1) for _ in range(u_num + 1)]
            for i in range(1, u_num + 1):
                for t in range(1, t_num + 1):
                    self.p[i][t] = m.getVarByName("p[%d,%d]" % (i, t)).x
            self.u = [[0] * (t_num + 1) for _ in range(u_num + 1)]
            for i in range(1, u_num + 1):
                for t in range(1, t_num + 1):
                    self.u[i][t] = m.getVarByName("u[%d,%d]" % (i, t)).x

        except gp.GurobiError as e:
            print('Error code ' + str(e))
        except AttributeError:
            print('Encountered an attribute error')

    def __init__(self, file_path, gurobi_params):
        """
        :param file_path: 求解样例文件的路径
        :param gurobi_params: 包含参数和参数值的字典
        """
        self.file_path = file_path
        self.gurobi_params = gurobi_params
        self.ObjVal = 0.0
        self.MIPGap = 0.0
        self.Runtime = 0.0
        self.p = []
        self.u = []
        self.solve(file_path)

    def write_in_file(self, file):
        """
        将求解结果写入文件
        :param file: Out/file
        :return:
        """
        with open("Out/" + file, "a+") as f:
            print("ObjVal: %f, Gap: %f, Runtime: %f" %
                  (self.ObjVal, self.MIPGap, self.Runtime), file=f)


s = MIQPSolver("../Data/5_std.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1})

# for i in range(1, 21):
#     solve("../Data/Based8Std/c%d_based_8_std.mod" % i)

# for i in range(12, 13):
#     solve("../Data/200Unit/200_0_%d_w.mod" % i)

# s.solve("")

# solve("../Data/8_std.mod")

# solve("../Data/10_0_2_w.mod")

# solve("../Data/200Unit/200_0_3_w.mod")
