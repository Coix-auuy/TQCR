f1 = open("../Out/QP/tqcr_p.out", 'r')
f2 = open("../Out/QP/p_tqcr_p.out", 'r')

for i in range(20):
    x = list(f1.readline().split())
    y = list(f2.readline().split())
    x = [float(i) for i in x]
    y = [float(i) for i in y]
    for j in range(0, len(x)):
        if abs(x[j] - y[j]) > 1e-4:
            print("yes")
            print(x[j])
            print(y[j])
            break

f1.close()
f2.close()
