class DataReader:

    def read_data(self, file_path):
        with open(file_path) as file:
            # 忽略第一行
            file.readline()
            # 读取第二行数据（运行时段）
            self.t_num = int(file.readline().split()[1])
            # 读取第三行数据（机组数量）
            self.u_num = int(file.readline().split()[1])
            # 忽略 4-10 行
            for _ in range(7):
                file.readline()
            # 读取第 11 行（负荷）
            self.load = list(file.readline().split())
            self.load = [float(_) for _ in self.load]
            # 忽略第 12 行
            file.readline()
            # 读取第 13 行（旋转保留）
            self.reserve = list(file.readline().split())
            self.reserve = [float(_) for _ in self.reserve]
            # 忽略第 14 行
            file.readline()
            # 读取机组参数
            u_params = []
            for _ in range(self.u_num):
                u_params.append(file.readline().split())
                temp = file.readline()
                self.r_up.append(float(temp.split()[1]))
                self.r_down.append(float(temp.split()[2]))
            # au + bp + cp^2
            self.alpha = [float(_[3]) for _ in u_params]
            self.beta = [float(_[2]) for _ in u_params]
            self.gama = [float(_[1]) for _ in u_params]
            # 最小出力
            self.min_output = [float(_[4]) for _ in u_params]
            # 最大出力
            self.max_output = [float(_[5]) for _ in u_params]
            # 初始运行时间
            self.t_init = [int(float(_[6])) for _ in u_params]
            # 初始出力
            self.p_init = [float(_[15]) for _ in u_params]
            # 最小开机时间
            self.t_on = [int(float(_[7])) for _ in u_params]
            # 最小停机时间
            self.t_off = [int(float(_[8])) for _ in u_params]
            # 热启动成本
            self.c_h = [float(_[13]) for _ in u_params]
            # 冷启动成本
            self.c_c = self.c_h
            if "std" in file_path:
                # 冷却时间
                self.t_c = [int(float(_[16])) for _ in u_params]
                if "5_std" not in file_path:
                    self.c_c = [2 * _ for _ in self.c_c]
            else:
                self.t_c = [1 for _ in range(self.u_num)]
            # 开机爬坡
            self.r_on = self.min_output
            # 停机爬坡
            self.r_off = self.min_output
            # u0
            self.u0 = [1 if _ > 0 else 0 for _ in self.p_init]
            if "8_std" in file_path:
                temp = [0.71, 0.65, 0.62, 0.6, 0.58, 0.58, 0.6, 0.64,
                        0.73, 0.8, 0.82, 0.83, 0.82, 0.8, 0.79, 0.79,
                        0.83, 0.91, 0.9, 0.88, 0.85, 0.84, 0.79, 0.74]

                self.load = [sum(self.max_output) * _ for _ in temp]
                self.reserve = [_ * 0.03 for _ in self.load]
        file.close()

    def __init__(self, file_path):

        # 运行时段
        self.t_num = 0
        # 机组数量
        self.u_num = 0
        # 每个时段的负荷
        self.load = []
        # 旋转保留
        self.reserve = []
        # 上爬坡
        self.r_up = []
        #  下爬坡
        self.r_down = []
        # 成本二次函数系数
        self.alpha = []
        self.beta = []
        self.gama = []
        # 发电功率下限
        self.min_output = []
        # 发电功率上限
        self.max_output = []
        # 初始状态持续了多久, (+) -> online, (-) -> offline
        self.t_init = []
        # 机组 0 时刻发电功率
        self.p_init = []
        # 最小开机时间
        self.t_on = []
        # 最小停机时间
        self.t_off = []
        # 热启动费用
        self.c_h = []
        # 冷启动费用
        self.c_c = []
        # 机组冷却时间
        self.t_c = []
        # 开机爬坡
        self.r_on = []
        # 停机爬坡
        self.r_off = []
        # 机组时刻 0 的状态
        self.u0 = []
        self.read_data(file_path)


s = DataReader("../Data/8_std.mod")
