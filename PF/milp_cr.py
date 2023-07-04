import gurobipy as gp


def solve(file_path):
    # Read data start
    with open(file_path) as file:
        # 忽略第一行
        file.readline()
        # 读取第二行数据（运行时段）
        period_num = int(file.readline().split()[1])
        # 读取第三行数据（机组数量）
        unit_num = int(file.readline().split()[1])
        # 忽略4~10行
        for i in range(7):
            file.readline()
        # 读取第11行数据（每个时段的负荷）
        load = list(file.readline().split())
        load = [float(i) for i in load]
        # 忽略第12行
        file.readline()
        # 读取第13行数据
        spinning_reserve = list(file.readline().split())  # 旋转备用
        spinning_reserve = [float(i) for i in spinning_reserve]
        # 忽略第14行
        file.readline()
        # 读取机组参数
        unit_parameter = []
        ramp_up = []  # 上爬坡
        ramp_down = []  # 下爬坡
        for i in range(unit_num):
            unit_parameter.append(file.readline().split())
            temp = file.readline()
            ramp_up.append(float(temp.split()[1]))
            ramp_down.append(float(temp.split()[2]))
        alpha = [float(i[3]) for i in unit_parameter]
        beta = [float(i[2]) for i in unit_parameter]
        gama = [float(i[1]) for i in unit_parameter]
        minimum_power_output = [float(i[4]) for i in unit_parameter]  # 发电功率下限
        maximum_power_output = [float(i[5]) for i in unit_parameter]  # 发电功率上限
        time_init = [float(i[6]) for i in unit_parameter]  # 到初始(0)状态持续了多久
        time_init = [int(i) for i in time_init]
        init_power_output = [float(i[15])
                             for i in unit_parameter]  # 机组初始0状态发电功率
        minimum_time_on = [float(i[7]) for i in unit_parameter]  # 最小开机时间
        minimum_time_on = [int(i) for i in minimum_time_on]
        minimum_time_off = [float(i[8]) for i in unit_parameter]  # 最小停机时间
        minimum_time_off = [int(i) for i in minimum_time_off]
        startup_cost_hot = [float(i[13]) for i in unit_parameter]  # 热启动费用
        startup_cost_cold = startup_cost_hot  # 冷启动费用，默认设为两倍的热启动费用
        if 'std' in file_path:
            hot_to_cold = [float(i[16]) for i in unit_parameter]  # 冷却时间
            hot_to_cold = [int(i) for i in hot_to_cold]
            if '5_std' not in file_path:
                startup_cost_cold = [2 * i for i in startup_cost_cold]  # 冷启动费用
        else:
            hot_to_cold = [1 for i in range(1, unit_num + 1)]
        startup_power_output = minimum_power_output  # 开机爬坡
        shutdown_power_output = minimum_power_output  # 停机爬坡
        # 机组时刻0的状态，根据0时刻机组出力得到
        u_i_0 = [1 if i > 0 else 0 for i in init_power_output]
        if '8_std' in file_path:
            temp = [0.71, 0.65, 0.62, 0.6, 0.58, 0.58, 0.6, 0.64,
                    0.73, 0.8, 0.82, 0.83, 0.82, 0.8, 0.79, 0.79,
                    0.83, 0.91, 0.9, 0.88, 0.85, 0.84, 0.79, 0.74]
            load = [sum(maximum_power_output) * i for i in temp]
            spinning_reserve = [i * 0.03 for i in load]
        file.close()
    # Read data end

    try:
        # Model start
        m = gp.Model("UC")

        # Create variables
        rc = gp.tuplelist([(i, t) for i in range(1, unit_num + 1)
                           for t in range(1, period_num + 1)])
        u = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="u")
        v = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="v")
        w = m.addVars(rc, lb=0.0, ub=1.0, vtype=gp.GRB.CONTINUOUS, name="w")
        p = m.addVars(rc,  vtype=gp.GRB.CONTINUOUS, name="p")
        sc = m.addVars(rc,  vtype=gp.GRB.CONTINUOUS, name="sc")
        z = m.addVars(rc, vtype=gp.GRB.CONTINUOUS, name="z")
        # 构造一些用于线性化的点
        l = 100
        c_p = [[] for _ in range(unit_num)]  # 构造出这些点
        for i in range(unit_num):
            for j in range(l + 1):
                c_p[i].append(minimum_power_output[i] + j *
                              (maximum_power_output[i] - minimum_power_output[i])/l)
        # Objective function (1)
        m.setObjective(gp.quicksum(
            z[i, t] + sc[i, t] for i in range(1, unit_num + 1) for t in
            range(1, period_num + 1)), gp.GRB.MINIMIZE)

        # Constraints
        # PF linearization constraints
        m.addConstrs((z[i, t] >= (2 * gama[i - 1]*c_p[i - 1][j] + beta[i - 1]) * p[i, t] + (alpha[i - 1] - gama[i - 1] * (c_p[i - 1][j] ** 2)) * u[i, t]
                     for i in range(1, unit_num+1) for t in range(1, period_num + 1) for j in range(1, l+1)), name="PF linearization constraints")

        # startup cost (2)(3)
        m.addConstrs((sc[i, t] >= v[i, t] * startup_cost_hot[i - 1] for i in range(1, unit_num + 1) for t in range(1, period_num + 1)),
                     name="startup cost constraints 1")
        m.addConstrs(
            (sc[i, j] >= startup_cost_cold[i - 1] * (
                v[i, j] - gp.quicksum(w[i, k] for k in range(max(j - minimum_time_off[i - 1] - hot_to_cold[i - 1], 1), j)) - (
                    1 if j - minimum_time_off[i - 1] - hot_to_cold[i - 1] <= 0 and max(0, -time_init[i - 1]) < abs(
                        j - minimum_time_off[i - 1] - hot_to_cold[i - 1]) + 1 else 0))
             for i in range(1, unit_num + 1) for j in range(1, period_num + 1)), name="startup cost constraints 2")

        # power balance constraints (4)
        m.addConstrs((gp.quicksum(p[i, t] for i in range(1, unit_num + 1)) == load[t - 1] for t in range(1, period_num + 1)),
                     name="power balance constraints")

        # system spinning reserve requirements (5)
        m.addConstrs((
            gp.quicksum(u[i, t] * maximum_power_output[i - 1] for i in range(1, unit_num + 1)) >= load[t - 1] +
            spinning_reserve[t - 1] for t in range(1, period_num + 1)), name="spinning reserve requirements")

        # unit generation limits (6)
        m.addConstrs((p[i, t] >= u[i, t] * minimum_power_output[i - 1] for i in range(1, unit_num + 1) for t in range(1, period_num + 1)),
                     name="unit generation limits 1")

        # unit generation limits (6)
        m.addConstrs((p[i, t] <= u[i, t] * maximum_power_output[i - 1] for i in range(1, unit_num + 1) for t in range(1, period_num + 1)),
                     name="unit generation limits 2")

        # ramp rate limits (7) (8)
        m.addConstrs(
            (p[i, t] - (p[i, t - 1] if t > 1 else init_power_output[i - 1]) <= (u[i, t - 1] if t > 1 else u_i_0[i - 1]) * ramp_up[i - 1] + v[i, t] *
             startup_power_output[i - 1] for i in range(1, unit_num + 1) for t in range(1, period_num + 1)), name="ramp rate limits 1")
        m.addConstrs(
            ((p[i, t - 1] if t > 1 else init_power_output[i - 1]) - p[i, t] <= u[i, t] * ramp_down[i - 1] + w[i, t] * shutdown_power_output[i - 1] for
             i in
             range(1, unit_num + 1) for t in range(1, period_num + 1)), name="ramp rate limits 2")

        # minimum up and down time constraints can be formulated (9)(10)
        W = []
        for i in range(1, unit_num + 1):
            W.append(max(
                0, min(period_num, u_i_0[i - 1] * (minimum_time_on[i - 1] - time_init[i - 1]))))
        L = []
        for i in range(1, unit_num + 1):
            L.append(max(0, min(
                period_num, (1 - u_i_0[i - 1]) * (minimum_time_off[i - 1] + time_init[i - 1]))))
        m.addConstrs(
            (gp.quicksum(v[i, j] for j in range(max(0, k - minimum_time_on[i - 1]) + 1, k + 1)) <= u[i, k] for i in range(1, unit_num + 1) for k
             in range(W[i - 1] + 1, period_num + 1)), name=" minimum on time constraints")
        m.addConstrs(
            (gp.quicksum(w[i, j] for j in range(max(0, k - minimum_time_off[i - 1]) + 1, k + 1)) <= 1 - u[i, k] for i in range(1, unit_num + 1)
             for k in range(L[i - 1] + 1, period_num + 1)), name="maximum off time constraints")

        #  state variables and logical constraints (11)
        m.addConstrs(
            (v[i, j] - w[i, j] == u[i, j] - (u[i, j - 1] if j > 1 else u_i_0[i - 1]) for i in range(1, unit_num + 1) for j in
             range(1, period_num + 1)),
            name="logical constraints")

        # initial status of units
        m.addConstrs(
            (u[i, j] == u_i_0[i - 1] for i in range(1, unit_num + 1)
             for j in range(1, W[i - 1] + L[i - 1] + 1)),
            name="initial status constraints")
        # Model end

        # Compute optimal cost
        m.update()
        m.optimize()
        objval = m.getAttr("ObjVal")
        f = open("Output/milp_cr.out", "a+")
        print(objval, file=f)
        f.close
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    except AttributeError:
        print('Encountered an attribute error')


for i in range(1, 21):
    solve("../Data/Based8Std/c%d_based_8_std.mod" % i)

# solve("../Data/8_std.mod")
