import gurobipy as gp
from common.data_reader import DataReader


class Solver:

    def solve(self):
        """
        调用求解器，子类重写
        :return:
        """
        pass

    def __init__(self, file_path, gurobi_params, org_or_cr, org_or_proj):
        """
        :param file_path: 求解样例文件的路径
        :param gurobi_params: 包含参数和参数值的字典
        :param org_or_cr: 0 - MIP or 1 - cr
        :param org_or_proj: 0 - MIP or 1 - P-MIP
        """
        self.file_path = file_path
        self.data = DataReader(self.file_path)
        self.gurobi_params = gurobi_params
        self.org_or_cr = org_or_cr
        self.org_or_proj = org_or_proj
        self.ObjVal = 0.0
        if self.org_or_cr == 0:
            self.MIPGap = 0.0
        self.Runtime = 0.0
        self.p = []
        self.u = []
        self.sc = []
        self.solve()

    def write_in_file(self, file):
        """
        将求解结果写入文件
        :param file: Out/file
        :return:
        """
        with open("Out/" + file, "a+") as f:
            if self.org_or_cr == 0:
                print("ObjVal: %f, Gap: %f, Runtime: %f" % (self.ObjVal, self.MIPGap, self.Runtime), file=f)
            else:
                print("ObjVal: %f, Runtime: %f" % (self.ObjVal, self.Runtime), file=f)


class QPSolver(Solver):
    def solve(self):
        """
        求解
        :return:
        """
        data = self.data
        t_num = data.t_num
        u_num = data.u_num
        load = data.load
        reserve = data.reserve
        r_up = data.r_up
        r_down = data.r_down
        alpha = data.alpha
        beta = data.beta
        gama = data.gama
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
        # 投影模型，相应变量做变换
        if self.org_or_proj == 1:
            alpha = [alpha[i] + beta[i] * min_output[i] + gama[i] * (min_output[i] ** 2) for i in range(u_num)]
            beta = [(max_output[i] - min_output[i]) * (beta[i] + 2 * gama[i] * min_output[i]) for i in range(u_num)]
            gama = [gama[i] * ((max_output[i] - min_output[i]) ** 2) for i in range(u_num)]
            r_up = [r_up[i] / (max_output[i] - min_output[i]) for i in range(u_num)]
            r_down = [r_down[i] / (max_output[i] - min_output[i]) for i in range(u_num)]
            p_init = [(p_init[i] - u0[i] * min_output[i]) / (max_output[i] - min_output[i]) for i in range(u_num)]
            r_on = [(r_on[i] - min_output[i]) / (max_output[i] - min_output[i]) for i in range(u_num)]
            r_off = [(r_off[i] - min_output[i]) / (max_output[i] - min_output[i]) for i in range(u_num)]

        try:
            # Model start
            m = gp.Model("QP")
            # Create variables
            rc = list([(i, t) for i in range(1, u_num + 1) for t in range(1, t_num + 1)])
            # org_or_cr == 1 : CR
            if self.org_or_cr == 1:
                u = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="u")
                v = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="v")
                w = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="w")
            else:
                u = m.addVars(rc, vtype=gp.GRB.BINARY, name="u")
                v = m.addVars(rc, vtype=gp.GRB.BINARY, name="v")
                w = m.addVars(rc, vtype=gp.GRB.BINARY, name="w")
            p = m.addVars(rc, lb=0, vtype=gp.GRB.CONTINUOUS, name="p")
            sc = m.addVars(rc, vtype=gp.GRB.CONTINUOUS, name="sc")

            # Objective function (1)
            m.setObjective(gp.quicksum(
                alpha[i - 1] * u[i, t] + beta[i - 1] * p[i, t] + gama[i - 1] * (p[i, t] ** 2) + sc[i, t] for i in range(1, u_num + 1) for t in
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

            # power balance constraints (4) or (23)
            if self.org_or_proj == 0:
                m.addConstrs((gp.quicksum(p[i, t] for i in range(1, u_num + 1)) == load[t - 1] for t in range(1, t_num + 1)),
                             name="power balance constraints")
            elif self.org_or_proj == 1:
                m.addConstrs((gp.quicksum(
                    p[i, t] * (max_output[i - 1] - min_output[i - 1]) + u[i, t] * min_output[i - 1] for i in range(1, u_num + 1)) == load[t - 1] for t
                              in range(1, t_num + 1)), name="power balance constraints")

            # system spinning reserve requirements (5)
            m.addConstrs((
                gp.quicksum(u[i, t] * max_output[i - 1] for i in range(1, u_num + 1)) >= load[t - 1] +
                reserve[t - 1] for t in range(1, t_num + 1)), name="spinning reserve requirements")

            # unit generation limits (6) or (24)
            if self.org_or_proj == 0:
                m.addConstrs((p[i, t] >= u[i, t] * min_output[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                             name="unit generation limits 1")
                m.addConstrs((p[i, t] <= u[i, t] * max_output[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                             name="unit generation limits 2")
            else:
                m.addConstrs((p[i, t] <= u[i, t] for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="unit generation limits")

            # ramp rate limits (7) (8) or (25) (26)
            m.addConstrs(
                (p[i, t] - (p[i, t - 1] if t > 1 else p_init[i - 1]) <= (u[i, t - 1] if t > 1 else u0[i - 1]) * r_up[i - 1] + v[i, t] *
                 r_on[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="ramp rate limits 1")
            m.addConstrs(
                ((p[i, t - 1] if t > 1 else p_init[i - 1]) - p[i, t] <= u[i, t] * r_down[i - 1] + w[i, t] * r_off[i - 1] for i in range(1, u_num + 1)
                 for t in range(1, t_num + 1)), name="ramp rate limits 2")

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
            if self.org_or_cr == 0:
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


# for i in range(1, 21):
#     s = QPSolver("../Data/Based8Std/c%d_based_8_std.mod" % i, {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 0, "Threads": 1}, 0)
#     del s

# for i in range(1, 13):
#     s = QPSolver("../Data/200Unit/200_0_%d_w.mod" % i, {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 0)
#     del s

# s = QPSolver("../Data/5_std.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1)

# s = QPSolver("../Data/8_std.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1, 0)

# s = QPSolver("../Data/10_0_2_w.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1, 1)

# s = QPSolver("../Data/200Unit/200_0_3_w.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1, 1)

# s = QPSolver("../Data/Based8Std/c13_based_8_std.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1, 1)

# print(s.ObjVal)

class SOCPSolver(Solver):
    def solve(self):
        """
        求解
        :return:
        """
        data = self.data
        t_num = data.t_num
        u_num = data.u_num
        load = data.load
        reserve = data.reserve
        r_up = data.r_up
        r_down = data.r_down
        alpha = data.alpha
        beta = data.beta
        gama = data.gama
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
        # 投影模型，相应变量做变换
        if self.org_or_proj == 1:
            alpha = [alpha[i] + beta[i] * min_output[i] + gama[i] * (min_output[i] ** 2) for i in range(u_num)]
            beta = [(max_output[i] - min_output[i]) * (beta[i] + 2 * gama[i] * min_output[i]) for i in range(u_num)]
            gama = [gama[i] * ((max_output[i] - min_output[i]) ** 2) for i in range(u_num)]
            r_up = [r_up[i] / (max_output[i] - min_output[i]) for i in range(u_num)]
            r_down = [r_down[i] / (max_output[i] - min_output[i]) for i in range(u_num)]
            p_init = [(p_init[i] - u0[i] * min_output[i]) / (max_output[i] - min_output[i]) for i in range(u_num)]
            r_on = [(r_on[i] - min_output[i]) / (max_output[i] - min_output[i]) for i in range(u_num)]
            r_off = [(r_off[i] - min_output[i]) / (max_output[i] - min_output[i]) for i in range(u_num)]

        try:
            # Model start
            m = gp.Model("SOCP")
            # Create variables
            rc = list([(i, t) for i in range(1, u_num + 1) for t in range(1, t_num + 1)])
            # org_or_cr == 1 : CR
            if self.org_or_cr == 1:
                u = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="u")
                v = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="v")
                w = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="w")
            else:
                u = m.addVars(rc, vtype=gp.GRB.BINARY, name="u")
                v = m.addVars(rc, vtype=gp.GRB.BINARY, name="v")
                w = m.addVars(rc, vtype=gp.GRB.BINARY, name="w")
            p = m.addVars(rc, vtype=gp.GRB.CONTINUOUS, name="p")
            sc = m.addVars(rc, vtype=gp.GRB.CONTINUOUS, name="sc")
            z = m.addVars(rc, vtype=gp.GRB.CONTINUOUS, name="z")

            # Objective function (1)
            m.setObjective(gp.quicksum(
                alpha[i - 1] * u[i, t] + beta[i - 1] * p[i, t] + gama[i - 1] * z[i, t] + sc[i, t] for i in range(1, u_num + 1) for t in
                range(1, t_num + 1)), gp.GRB.MINIMIZE)

            # Constraints
            # p^2 <= uz
            m.addConstrs((p[i, t] ** 2 <= u[i, t] * z[i, t] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                         name="quadratic constraints 1")
            m.addConstrs((z[i, t] >= 0 for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="quadratic constraints 2")
            # m.addConstrs((u[i, t] * z[i, t] <= max_output[i - 1] ** 2 for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
            #              name="quadratic constraints 3")
            # startup cost (2)(3)
            m.addConstrs((sc[i, t] >= v[i, t] * c_h[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                         name="startup cost constraints 1")
            m.addConstrs(
                (sc[i, t] >= c_c[i - 1] * (v[i, t] - gp.quicksum(w[i, k] for k in range(max(t - t_off[i - 1] - t_c[i - 1], 1), t)) - (
                    1 if t - t_off[i - 1] - t_c[i - 1] <= 0 and max(0, -t_init[i - 1]) < abs(
                        t - t_off[i - 1] - t_c[i - 1] - 1) + 1 else 0))
                 for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="startup cost constraints 2")

            # power balance constraints (4) or (23)
            if self.org_or_proj == 0:
                m.addConstrs((gp.quicksum(p[i, t] for i in range(1, u_num + 1)) == load[t - 1] for t in range(1, t_num + 1)),
                             name="power balance constraints")
            elif self.org_or_proj == 1:
                m.addConstrs((gp.quicksum(
                    p[i, t] * (max_output[i - 1] - min_output[i - 1]) + u[i, t] * min_output[i - 1] for i in range(1, u_num + 1)) == load[t - 1] for t
                              in range(1, t_num + 1)), name="power balance constraints")

            # system spinning reserve requirements (5)
            m.addConstrs((
                gp.quicksum(u[i, t] * max_output[i - 1] for i in range(1, u_num + 1)) >= load[t - 1] +
                reserve[t - 1] for t in range(1, t_num + 1)), name="spinning reserve requirements")

            # unit generation limits (6) or (24)
            if self.org_or_proj == 0:
                m.addConstrs((p[i, t] >= u[i, t] * min_output[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                             name="unit generation limits 1")
                m.addConstrs((p[i, t] <= u[i, t] * max_output[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                             name="unit generation limits 2")
            else:
                m.addConstrs((p[i, t] <= u[i, t] for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="unit generation limits 1")
                m.addConstrs((p[i, t] >= 0 for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="unit generation limits 2")

            # ramp rate limits (7) (8) or (25) (26)
            m.addConstrs(
                (p[i, t] - (p[i, t - 1] if t > 1 else p_init[i - 1]) <= (u[i, t - 1] if t > 1 else u0[i - 1]) * r_up[i - 1] + v[i, t] *
                 r_on[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="ramp rate limits 1")
            m.addConstrs(
                ((p[i, t - 1] if t > 1 else p_init[i - 1]) - p[i, t] <= u[i, t] * r_down[i - 1] + w[i, t] * r_off[i - 1] for i in range(1, u_num + 1)
                 for t in range(1, t_num + 1)), name="ramp rate limits 2")

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
            if self.org_or_cr == 0:
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
            self.sc = [[0] * (t_num + 1) for _ in range(u_num + 1)]
            for i in range(1, u_num + 1):
                for t in range(1, t_num + 1):
                    self.sc[i][t] = m.getVarByName("sc[%d,%d]" % (i, t)).x

        except gp.GurobiError as e:
            print('Error code ' + str(e))
        except AttributeError:
            print('Encountered an attribute error')


# for i in range(1, 21):
#     s = SOCPSolver("../Data/Based8Std/c%d_based_8_std.mod" % i, {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 0, "Threads": 1}, 0)
#     del s

# for i in range(1, 13):
#     s = SOCPSolver("../Data/200Unit/200_0_%d_w.mod" % i, {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 0)
#     del s

# s = SOCPSolver("../Data/5_std.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1, 0)

# s = SOCPSolver("../Data/8_std.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 0, 1)


# s = SOCPSolver("../Data/10_0_2_w.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1}, 1, 0)

# s = SOCPSolver("../Data/200Unit/200_0_3_w.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1, 0)

# s = SOCPSolver("../Data/Based8Std/c13_based_8_std.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1, 1)

# print(s.ObjVal)


class QCRSolver(Solver):
    def solve(self):
        """
        求解
        :return:
        """
        # 拿到样例数据
        data = self.data
        t_num = data.t_num
        u_num = data.u_num
        load = data.load
        reserve = data.reserve
        r_up = data.r_up
        r_down = data.r_down
        alpha = data.alpha
        beta = data.beta
        gama = data.gama
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
        # 投影模型，相应变量做变换
        if self.org_or_proj == 1:
            alpha = [alpha[i] + beta[i] * min_output[i] + gama[i] * (min_output[i] ** 2) for i in range(u_num)]
            beta = [(max_output[i] - min_output[i]) * (beta[i] + 2 * gama[i] * min_output[i]) for i in range(u_num)]
            gama = [gama[i] * ((max_output[i] - min_output[i]) ** 2) for i in range(u_num)]
            r_up = [r_up[i] / (max_output[i] - min_output[i]) for i in range(u_num)]
            r_down = [r_down[i] / (max_output[i] - min_output[i]) for i in range(u_num)]
            p_init = [(p_init[i] - u0[i] * min_output[i]) / (max_output[i] - min_output[i]) for i in range(u_num)]
            r_on = [(r_on[i] - min_output[i]) / (max_output[i] - min_output[i]) for i in range(u_num)]
            r_off = [(r_off[i] - min_output[i]) / (max_output[i] - min_output[i]) for i in range(u_num)]
        # 计算 QCR 附加项的系数
        socp_cr = SOCPSolver(self.file_path, {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": self.gurobi_params["OutputFlag"], "Threads": 1}, 1,
                             self.org_or_proj)
        upsilon = [[0] * (t_num + 1) for _ in range(u_num + 1)]
        for i in range(1, u_num + 1):
            for t in range(1, t_num + 1):
                if socp_cr.u[i][t] == 0:
                    upsilon[i][t] = 0
                else:
                    upsilon[i][t] = -2 * gama[i - 1] * socp_cr.p[i][t] / socp_cr.u[i][t]

        tau = [[0] * (t_num + 1) for _ in range(u_num + 1)]
        for i in range(1, u_num + 1):
            for t in range(1, t_num + 1):
                if socp_cr.u[i][t] == 0:
                    tau[i][t] = 0
                else:
                    tau[i][t] = gama[i - 1] * ((socp_cr.p[i][t] / socp_cr.u[i][t]) ** 2)
        # 解决二次项目标系数范围太大的问题
        for i in range(1, u_num + 1):
            for t in range(1, t_num + 1):
                if abs(upsilon[i][t]) <= 1e-5:
                    tau[i][t] = 0
                    upsilon[i][t] = 0
                if tau[i][t] <= 1e-5:
                    tau[i][t] = 0
                    upsilon[i][t] = 0

        try:
            # Model start
            m = gp.Model("SOCP")
            # Create variables
            rc = list([(i, t) for i in range(1, u_num + 1) for t in range(1, t_num + 1)])
            # org_or_cr == 1 : CR
            if self.org_or_cr == 1:
                u = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="u")
                v = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="v")
                w = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="w")
            else:
                u = m.addVars(rc, vtype=gp.GRB.BINARY, name="u")
                v = m.addVars(rc, vtype=gp.GRB.BINARY, name="v")
                w = m.addVars(rc, vtype=gp.GRB.BINARY, name="w")
            p = m.addVars(rc, vtype=gp.GRB.CONTINUOUS, name="p")
            sc = m.addVars(rc, vtype=gp.GRB.CONTINUOUS, name="sc")

            # Objective function (1)
            m.setObjective(gp.quicksum(
                alpha[i - 1] * u[i, t] + beta[i - 1] * p[i, t] + gama[i - 1] * (p[i, t] ** 2) + sc[i, t] + tau[i][t] * (u[i, t] ** 2 - u[i, t]) +
                upsilon[i][t] * (u[i, t] * p[i, t] - p[i, t]) for i in range(1, u_num + 1) for t in
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

            # power balance constraints (4) or (23)
            if self.org_or_proj == 0:
                m.addConstrs((gp.quicksum(p[i, t] for i in range(1, u_num + 1)) == load[t - 1] for t in range(1, t_num + 1)),
                             name="power balance constraints")
            elif self.org_or_proj == 1:
                m.addConstrs((gp.quicksum(
                    p[i, t] * (max_output[i - 1] - min_output[i - 1]) + u[i, t] * min_output[i - 1] for i in range(1, u_num + 1)) == load[t - 1] for t
                              in range(1, t_num + 1)), name="power balance constraints")

            # system spinning reserve requirements (5)
            m.addConstrs((
                gp.quicksum(u[i, t] * max_output[i - 1] for i in range(1, u_num + 1)) >= load[t - 1] +
                reserve[t - 1] for t in range(1, t_num + 1)), name="spinning reserve requirements")

            # unit generation limits (6) or (24)
            if self.org_or_proj == 0:
                m.addConstrs((p[i, t] >= u[i, t] * min_output[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                             name="unit generation limits 1")
                m.addConstrs((p[i, t] <= u[i, t] * max_output[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)),
                             name="unit generation limits 2")
            else:
                m.addConstrs((p[i, t] <= u[i, t] for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="unit generation limits 1")
                m.addConstrs((p[i, t] >= 0 for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="unit generation limits 2")

            # ramp rate limits (7) (8) or (25) (26)
            m.addConstrs(
                (p[i, t] - (p[i, t - 1] if t > 1 else p_init[i - 1]) <= (u[i, t - 1] if t > 1 else u0[i - 1]) * r_up[i - 1] + v[i, t] *
                 r_on[i - 1] for i in range(1, u_num + 1) for t in range(1, t_num + 1)), name="ramp rate limits 1")
            m.addConstrs(
                ((p[i, t - 1] if t > 1 else p_init[i - 1]) - p[i, t] <= u[i, t] * r_down[i - 1] + w[i, t] * r_off[i - 1] for i in range(1, u_num + 1)
                 for t in range(1, t_num + 1)), name="ramp rate limits 2")

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
            if self.org_or_cr == 0:
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

# for i in range(1, 21):
#     s = QCRSolver("../Data/Based8Std/c%d_based_8_std.mod" % i, {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 0, "Threads": 1}, 0)
#     del s

# for i in range(1, 13):
#     s = QCRSolver("../Data/200Unit/200_0_%d_w.mod" % i, {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 0)
#     del s

# s = QCRSolver("../Data/5_std.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1)

# s = QCRSolver("../Data/8_std.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 0, 1)

# s = QCRSolver("../Data/10_0_2_w.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1, 0)

# s = QCRSolver("../Data/200Unit/200_0_3_w.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1, 0)

# s = QCRSolver("../Data/Based8Std/c13_based_8_std.mod", {"MIPGap": 0, "TimeLimit": 900, "OutputFlag": 1, "Threads": 1}, 1, 0)

# print(s.ObjVal)
