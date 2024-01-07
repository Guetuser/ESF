import sys
from gurobipy import *
import numpy as np

np.set_printoptions(threshold=np.inf)


# S_box = [12, 5, 6, 10, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]   # present的 S盒
# S_box = [1, 10, 4, 12, 6, 15, 3, 9, 2, 13, 11, 7, 5, 0, 8, 14]  # GIFT 的 S盒
#S_box = [1, 2, 4, 13, 6, 15, 11, 8, 10, 5, 14, 3, 9, 12, 7, 0]  # PICO 的 S盒
#S_box = [12, 0, 15, 10, 2, 11, 9, 5, 8, 3, 13, 7, 1, 14, 6, 4]    # TWINE 的Sbox

# 文本读取sagemath生成的H表达式
def duqu_H(wenjian):
    temp = []
    with open(wenjian, 'r') as f:
        content = [line.rstrip('\n') for line in f]
        for line in content:
            line = line[:-4].replace('x +', '').replace(' (', '').replace(') ', '').replace(', ', ' ').replace(' ', ',')
            line = line[:-1].split(',')
            temp.append(line)
    return temp


# str-->int
def str2int(list):
    for num in list:
        for i in range(len(num)):
            num[i] = int(num[i])
    return list


# 计算结果
def calculate(list1, list2):
    H_list = np.array(list1)
    Impossible_list = np.array(list2)
    jieguo = H_list[:-1] * Impossible_list
    jieguo = sum(jieguo) + H_list[-1]
    # print(jieguo)
    return jieguo


# 贪心算法
def Greedy(H, X):
    # print(len(H), len(X))  #  建表
    biao = np.zeros(shape=(len(H), len(X)))
    for i in range(len(H)):
        for j in range(len(X)):
            if calculate(H[i], X[j]) < 0:
                biao[i][j] = 1
    # print(biao)
    # print(np.sum(biao, axis=0))
    one_num = []
    one_index = []
    for i in range(len(X)):
        for j in range(len(H)):
            one_num.append(int(biao[j][i]))
            if biao[j][i] == 1:
                one_index.append(j)
            else:
                one_index.append('None')
    one_index = np.array(one_index).reshape(len(X), len(H))
    one_index = one_index.tolist()

    for i in one_index:
        while 'None' in i:
            i.remove('None')
    # print(len(one_index[0]),one_index[0])

    one_num = np.array(one_num).reshape(len(X), len(H))
    # print(sum(one_num[0]), sum(one_num[3]), sum(one_num[-3]), sum(one_num[-1]))
    return one_index


# gurobi优化器
def gubi(file):
    model = read(file)

    if model.isMIP == 0:
        print('Model is not a MIP')
        exit(0)

    model.optimize()

    if model.status == GRB.Status.OPTIMAL:
        print('Optimal objective: %g' % model.objVal)
    elif model.status == GRB.Status.INF_OR_UNBD:
        print('Model is infeasible or unbounded')
        exit(0)
    elif model.status == GRB.Status.INFEASIBLE:
        print('Model is infeasible')
        exit(0)

    elif model.status == GRB.Status.UNBOUNDED:
        print('Model is unbounded')
        exit(0)
    else:
        print('Optimization ended with status %d' % model.status)
        exit(0)

    # Iterate over the solutions and compute the objectives
    model.Params.outputFlag = 0
    print('')
    for k in range(model.solCount):
        model.Params.solutionNumber = k
        objn = 0
        for v in model.getVars():
            objn += v.obj * v.xn
        print('Solution %d has objective %g' % (k, objn))
    print('')

    print('---------------------------------------')
    for v in model.getVars():
        if v.x ==1:
            print(v.varName, v.x)
    print('---------------------------------------')
    print('Obj:', model.ObjVal)

    print('---------------------------------------')
    model.Params.outputFlag = 1

    fixed = model.fixed()
    fixed.Params.presolve = 0
    fixed.optimize()

    if fixed.status != GRB.Status.OPTIMAL:
        print("Error: fixed model isn't optimal")
        exit(1)

    diff = model.objVal - fixed.objVal

    if abs(diff) > 1e-6 * (1.0 + abs(model.objVal)):
        print('Error: objective values are different')
        exit(1)
    result = []
    # Print values of nonzero variables
    for v in fixed.getVars():
        if v.x != 0:
            # print('%s %g' % (v.varName, v.x))
            # print(v.varName)
            result.append(v.varName.replace('z', ''))

    return result

#  第二次gurobi，找不可能差分路径
def gubi_2(file):
    final = []
    model = read(file)
    if model.isMIP == 0:
        print('Model is not a MIP')
        exit(0)

    model.optimize()

    if model.status == GRB.Status.OPTIMAL:
        print('Optimal objective: %g' % model.objVal)
    elif model.status == GRB.Status.INF_OR_UNBD:
        print('Model is infeasible or unbounded')
        exit(0)
    elif model.status == GRB.Status.INFEASIBLE:
        print('######################################################################')
        # print(i, j)
        # print('######################################################################')
        # print('Model is infeasible')
        # final.append(i)
        # final.append(j)

        # exit(0)

    elif model.status == GRB.Status.UNBOUNDED:
        print('Model is unbounded')
        exit(0)
    else:
        print('Optimization ended with status %d' % model.status)
        exit(0)

        # Iterate over the solutions and compute the objectives
        model.Params.outputFlag = 0
        print('')
        for k in range(model.solCount):
            model.Params.solutionNumber = k
            objn = 0
            for v in model.getVars():
                objn += v.obj * v.xn
            print('Solution %d has objective %g' % (k, objn))
        print('')
        model.Params.outputFlag = 1

        fixed = model.fixed()
        fixed.Params.presolve = 0
        fixed.optimize()

        if fixed.status != GRB.Status.OPTIMAL:
            print("Error: fixed model isn't optimal")
            exit(1)

        diff = model.objVal - fixed.objVal

        if abs(diff) > 1e-6 * (1.0 + abs(model.objVal)):
            print('Error: objective values are different')
            exit(1)
        result = []
        # Print values of nonzero variables
        for v in fixed.getVars():
            if v.x != 0:
                # print('%s %g' % (v.varName, v.x))
                # print(v.varName)
                result.append(v.varName.replace('z', ''))
        # return result
    return model.status


#
def Xor(i,input1,input2,output,dummy):
    buf = ''
    for j in range(32):
        buf = buf + str(input1) + str(i) + "_" + str((j+25)%32+32) + " + " + str(input2) + str(i) + "_" + str(j) + \
                " + " + str(output) + str(i+1) + "_" + str(j) + " - " + "2 "+ str(dummy) + str(i) + "_" + str(j) + " = 0 " + "\n"

    return buf
# #     # opOuter.write(buf)
#
def zhiwei(i,input,output):
    buf = ''
    for j in range(32):
        buf = buf + str(input) + str(i) + "_" + str(j) + " - " + str(output) + str(i+1) + "_" + str(j+32) + " = 0" + "\n"
    return buf
#
#     # opOuter.write(buf)




def PrintOuter(first, last):
    opOuter = open("Outer.lp", 'w+')
    # opOuter.write("Minimize\n")
    # buf = ''
    # for i in range(0, ROUND):
    #     for j in range(0, 8):
    #         buf = buf + "a" + str(i) + "_" + str(j)
    #         if i != ROUND - 1 or j != 7:
    #             buf = buf + " + "
    # opOuter.write(buf)
    # opOuter.write('\n')
#
    opOuter.write("Subject to\n")
    buf = ''

    buf = ''
    for j in range(0, 64):
        buf = buf + "x" + str(0) + "_" + str(j) + " = " + str(first[0][j]) + "\n"
    opOuter.write(buf)

    buf = ''
    for j in range(0, 64):
        buf = buf + "x" + str(ROUND) + "_" + str(j) + " = " + str(last[0][j]) + "\n"
    opOuter.write(buf)

    for i in range(0, ROUND):
        buf = ''
        for j in range(0, 8):
            buf = ''
            for k in range(0, 4):
                buf = buf + "x" + str(i) + "_" + str(4 * j + k)
                if k != 3:
                    buf = buf + " + "
            buf = buf + " - a" + str(i) + "_" + str(j) + " >= 0\n"
            for k in range(0, 4):
                buf = buf + "x" + str(i) + "_" + str(4 * j + k) + " - a" + str(i) + "_" + str(j) + " <= 0\n"

            if j == 0:
                for k in range(0, 21):
                    for l in range(0, 9):
                        if conv0[9 * k + l] > 0:
                            if l <= 3:
                                buf = buf + " + " + str(conv0[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l )
                            if 4 <= l and l <= 7:
                                buf = buf + " + " + str(conv0[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= -" + str(conv0[9 * k + l]) + "\n"
                        if conv0[9 * k + l] < 0:
                            if l <= 3:
                                buf = buf + " - " + str(-conv0[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l)
                            if 4 <= l and l <= 7:
                                buf = buf + " - " + str(-conv0[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= " + str(-conv0[9 * k + l]) + "\n"
                        if conv0[9 * k + l] == 0:
                            if l == 8:
                                buf = buf + " >= " + str(conv0[9 * k + l]) + "\n"
            if j == 1:
                for k in range(0, 21):
                    for l in range(0, 9):
                        if conv1[9 * k + l] > 0:
                            if l <= 3:
                                buf = buf + " + " + str(conv1[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l )
                            if 4 <= l and l <= 7:
                                buf = buf + " + " + str(conv1[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= -" + str(conv1[9 * k + l]) + "\n"
                        if conv1[9 * k + l] < 0:
                            if l <= 3:
                                buf = buf + " - " + str(-conv1[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l)
                            if 4 <= l and l <= 7:
                                buf = buf + " - " + str(-conv1[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= " + str(-conv1[9 * k + l]) + "\n"
                        if conv1[9 * k + l] == 0:
                            if l == 8:
                                buf = buf + " >= " + str(conv1[9 * k + l]) + "\n"
            if j == 2:
                for k in range(0, 21):
                    for l in range(0, 9):
                        if conv2[9 * k + l] > 0:
                            if l <= 3:
                                buf = buf + " + " + str(conv2[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l )
                            if 4 <= l and l <= 7:
                                buf = buf + " + " + str(conv2[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= -" + str(conv2[9 * k + l]) + "\n"
                        if conv2[9 * k + l] < 0:
                            if l <= 3:
                                buf = buf + " - " + str(-conv2[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l)
                            if 4 <= l and l <= 7:
                                buf = buf + " - " + str(-conv2[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= " + str(-conv2[9 * k + l]) + "\n"
                        if conv2[9 * k + l] == 0:
                            if l == 8:
                                buf = buf + " >= " + str(conv2[9 * k + l]) + "\n"
            if j == 3:
                for k in range(0, 27):
                    for l in range(0, 9):
                        if conv3[9 * k + l] > 0:
                            if l <= 3:
                                buf = buf + " + " + str(conv3[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l )
                            if 4 <= l and l <= 7:
                                buf = buf + " + " + str(conv3[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= -" + str(conv3[9 * k + l]) + "\n"
                        if conv3[9 * k + l] < 0:
                            if l <= 3:
                                buf = buf + " - " + str(-conv3[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l)
                            if 4 <= l and l <= 7:
                                buf = buf + " - " + str(-conv3[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= " + str(-conv3[9 * k + l]) + "\n"
                        if conv3[9 * k + l] == 0:
                            if l == 8:
                                buf = buf + " >= " + str(conv3[9 * k + l]) + "\n"
            if j == 4:
                for k in range(0, 23):
                    for l in range(0, 9):
                        if conv4[9 * k + l] > 0:
                            if l <= 3:
                                buf = buf + " + " + str(conv4[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l)
                            if 4 <= l and l <= 7:
                                buf = buf + " + " + str(conv4[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= -" + str(conv4[9 * k + l]) + "\n"
                        if conv4[9 * k + l] < 0:
                            if l <= 3:
                                buf = buf + " - " + str(-conv4[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l)
                            if 4 <= l and l <= 7:
                                buf = buf + " - " + str(-conv4[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= " + str(-conv4[9 * k + l]) + "\n"
                        if conv4[9 * k + l] == 0:
                            if l == 8:
                                buf = buf + " >= " + str(conv4[9 * k + l]) + "\n"
            if j == 5:
                for k in range(0, 23):
                    for l in range(0, 9):
                        if conv5[9 * k + l] > 0:
                            if l <= 3:
                                buf = buf + " + " + str(conv5[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l )
                            if 4 <= l and l <= 7:
                                buf = buf + " + " + str(conv5[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= -" + str(conv5[9 * k + l]) + "\n"
                        if conv5[9 * k + l] < 0:
                            if l <= 3:
                                buf = buf + " - " + str(-conv5[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l)
                            if 4 <= l and l <= 7:
                                buf = buf + " - " + str(-conv5[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= " + str(-conv5[9 * k + l]) + "\n"
                        if conv5[9 * k + l] == 0:
                            if l == 8:
                                buf = buf + " >= " + str(conv5[9 * k + l]) + "\n"
            if j == 6:
                for k in range(0, 21):
                    for l in range(0, 9):
                        if conv6[9 * k + l] > 0:
                            if l <= 3:
                                buf = buf + " + " + str(conv6[9 * k + l]) + " x" + str(i) + "_" + str(
                                    4 * j + 3 - l )
                            if 4 <= l and l <= 7:
                                buf = buf + " + " + str(conv6[9 * k + l]) + " m" + str(i) + "_" + str(
                                    P[4 * j + 7 - l ])
                            if l == 8:
                                buf = buf + " >= -" + str(conv6[9 * k + l]) + "\n"
                        if conv6[9 * k + l] < 0:
                            if l <= 3:
                                buf = buf + " - " + str(-conv6[9 * k + l]) + " x" + str(i) + "_" + str(
                                    4 * j + 3 - l )
                            if 4 <= l and l <= 7:
                                buf = buf + " - " + str(-conv6[9 * k + l]) + " m" + str(i) + "_" + str(
                                    P[4 * j + 7 - l ])
                            if l == 8:
                                buf = buf + " >= " + str(-conv6[9 * k + l]) + "\n"
                        if conv6[9 * k + l] == 0:
                            if l == 8:
                                buf = buf + " >= " + str(conv6[9 * k + l]) + "\n"
            if j == 7:
                for k in range(0, 27):
                    for l in range(0, 9):
                        if conv7[9 * k + l] > 0:
                            if l <= 3:
                                buf = buf + " + " + str(conv7[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l )
                            if 4 <= l and l <= 7:
                                buf = buf + " + " + str(conv7[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= -" + str(conv7[9 * k + l]) + "\n"
                        if conv7[9 * k + l] < 0:
                            if l <= 3:
                                buf = buf + " - " + str(-conv7[9 * k + l]) + " x" + str(i) + "_" + str(4 * j + 3 - l)
                            if 4 <= l and l <= 7:
                                buf = buf + " - " + str(-conv7[9 * k + l]) + " m" + str(i) + "_" + str(P[4 * j + 7 - l])
                            if l == 8:
                                buf = buf + " >= " + str(-conv7[9 * k + l]) + "\n"
                        if conv7[9 * k + l] == 0:
                            if l == 8:
                                buf = buf + " >= " + str(conv7[9 * k + l]) + "\n"
            opOuter.write(buf)
        opOuter.write(Xor(i, "x", "m", "x", "d"))
        opOuter.write(zhiwei(i, "x", "x"))

    opOuter.write("Binary\n")
    buf = ''
    for i in range(0, ROUND):
        buf = ''
        for j in range(0, 8):
            buf = buf + "a" + str(i) + "_" + str(j) + "\n"
        opOuter.write(buf)
    for i in range(0, ROUND + 1):
        buf = ''
        for j in range(0, 64):
            buf = buf + "x" + str(i) + "_" + str(j) + "\n"
        opOuter.write(buf)

    for i in range(0, ROUND):
        buf = ''
        for j in range(0, 32):
            buf = buf + "d" + str(i) + "_" + str(j) + "\n"
        opOuter.write(buf)

    for i in range(0, ROUND):
        buf = ''
        for j in range(32):
            buf = buf + "m" + str(i) + "_" + str(j) + "\n"
        opOuter.write(buf)
#
#
#
    opOuter.close()



if __name__ == '__main__':
    file = "greedy_FSE.lp"
    file_zuizhong = "Outer.lp"
    output_file = "imp_chanfen_FSE_9R.txt"
    # # P置换
    # # P = [5,0,1,4,7,12,3,8,13,6,9,2,15,10,11,14]   没有bit置换
    P = (0, 8, 16, 24, 1, 9, 17, 25, 2, 10, 18, 26, 3, 11, 19, 27, 4, 12, 20, 28, 5, 13, 21, 29, 6, 14, 22, 30, 7, 15, 23, 31)


    ROUND = 9 #
    # # act = (1, 2, 3, 5, 7, 10, 13, 16, 18, 20, 22, 24, 26)
    # act = (0,1,1,2,3,4,6,8)
    BanListlen = 0

    # jieguo = gubi(file)
    # H = H_representation1()
    # print(H)

    conv0 = [0, 0, -1, -1, -1, 0, -1, 1, 3, -2, -1, 0, -2, -2, 2, -1, -1, 7, -1, 0, -1, 0, -1, 0, 1, -1, 3, 2, -1, -3, -2, -1, -2, -2, 1, 8, -2, 1, 0, -2, 3, 1, 4, 4, 0, 2, 1, 1, 2, 0, 1, -1, -1, 0, 2, -2, 2, 0, 1, -1, -2, -1, 4, 2, -1, 2, -1, 2, 3, 2, -1, 0, 3, 4, 1, -2, -2, 4, 1, 3, 0, -1, 1, 1, -1, 0, 0, -1, -1, 3, 2, 1, 1, -1, -2, -1, 1, -2, 4, 3, -1, -1, 3, -1, 2, 2, 2, 0, -1, -1, 2, 2, 2, 3, -1, 2, 0, 1, 2, -2, 1, 3, -1, 1, 1, 0, -2, 4, 1, 3, -2, 4, 3, 1, 0, -1, -2, 1, 1, 2, -2, 1, -2, 5, 1, -1, -1, 1, -1, 0, 0, 0, 2, 3, 1, 2, 3, 2, -4, 1, 1, 0, -2, 1, 2, 3, -4, -4, -4, -1, 11, -3, -3, 1, -2, -2, -4, 1, 3, 10, -1, -1, -2, -1, 2, 0, -2, -2, 7]
    conv1 = [-1, -2, -2, 1, -2, -1, 2, -3, 8, -1, -2, 1, -2, 2, -1, -2, -3, 8, -1, 2, -1, -1, -2, 0, -2, -1, 6, 2, -4, 1, 1, 3, 1, 3, 2, 0, 1, 0, -1, -1, -1, 1, -1, 0, 3, 2, 3, -1, 2, 2, -1, -1, 2, 0, 2, 3, 2, -1, -1, -1, 2, 2, 0, -1, -2, -2, 1, 1, 2, -2, -1, 6, 3, 1, 4, 4, -2, 1, -2, 0, 0, 1, -1, 3, -1, -2, -3, -3, 3, 7, -1, 2, 2, 2, 3, -1, 3, -1, 0, -1, -2, 1, -2, -2, 2, 1, -1, 6, -2, 4, 3, 1, 3, 4, -2, 1, 0, -2, 4, 1, 3, -2, 4, 3, 1, 0, 1, -1, -1, 3, -3, -3, -2, 3, 7, -1, 0, 0, 0, 1, -1, 1, -1, 2, 1, 1, -1, -1, 2, 2, 2, 0, 0, -2, -4, 1, 3, -2, -3, -3, 1, 10, -1, -2, -4, -4, 2, -1, 2, 3, 8, -4, -4, -1, -4, -2, 1, 3, 2, 11, 1, 1, -2, -2, -1, -1, -1, -2, 7]
    conv2 = [-2, 1, -2, 2, -1, 1, -1, -3, 6, -1, -1, 1, -1, -1, -1, -1, -2, 6, -1, 2, -2, -2, -1, 1, -2, 1, 6, -2, -1, -1, 5, 4, 4, 4, -2, 0, 1, -1, 1, -1, -2, 1, 1, -2, 4, 1, 1, 1, 4, -2, -2, 2, 3, 0, 3, 2, 3, 1, 3, 2, -4, 1, 0, 0, -1, -1, 1, -1, 1, -1, 0, 3, 4, 2, 3, 4, 1, -1, -1, -2, 0, -1, -2, 2, -1, -2, 2, -1, 1, 5, -2, 1, -2, -1, 1, -3, 2, -3, 8, 1, 1, -1, -2, 1, 1, -2, -2, 5, 2, -2, 1, -1, -2, -1, -2, 1, 6, 1, 1, -1, -1, -2, -3, 3, -3, 7, 1, 2, 1, -3, 2, 1, 2, 2, 0, -1, 2, 3, 1, 1, -1, -1, 3, 0, 3, 2, -1, 3, 1, -1, -1, 2, 0, 1, -1, 1, 1, 1, -1, -1, -2, 3, -1, 1, 2, -1, 2, 0, 1, 1, 0, 1, -2, -3, -3, 1, 4, 3, 2, 4, -1, -2, -2, -1, 1, -2, -1, 2, 7]
    conv3 = [0, 1, -2, -2, 1, -2, -2, 1, 6, -2, -1, 2, -1, -2, -1, 1, 2, 5, 0, 1, -1, -1, -1, 1, -1, -1, 4, -1, -2, -2, 1, 1, -1, 1, 2, 4, -2, -1, -1, 1, 2, 1, -2, -2, 6, -3, -2, 3, 1, -2, 1, -1, -3, 8, 1, 2, 0, -1, -2, -1, -2, 2, 4, 1, 1, 2, -3, 2, 1, 2, 2, 0, 3, 2, 3, 3, -1, -1, -1, 0, 0, -1, 1, 2, -2, -2, 1, 2, -2, 5, -2, 2, -1, 1, -2, -2, -1, -1, 7, 2, 1, -2, 2, 2, 1, 0, 1, 0, -1, 3, 2, 3, -1, 3, -1, 2, 0, -1, 2, 2, 2, 3, -1, 2, -1, 0, 2, -1, 1, -2, 1, -2, -1, -1, 5, 1, -3, -1, -3, -2, -2, -1, -3, 12, -2, 1, -2, 2, 1, 4, 4, 3, 0, 1, -2, -1, 2, -2, 1, 1, -2, 5, 3, 1, -3, -1, -2, -3, 3, -1, 7, 1, -2, -2, 1, -1, 1, -2, 1, 5, 2, -1, 2, 2, 3, -1, 3, -1, 0, -1, 2, -1, -2, 1, -1, -2, -1, 6, -1, -2, 1, 1, 1, 2, -2, -1, 4, 1, -3, 1, 1, 2, 2, 1, 2, 0, -1, -2, -2, -2, 1, 1, 1, -1, 6, 1, -1, -2, -2, 4, 5, 5, 3, 0, -2, -1, 2, -2, -1, -1, -1, 1, 6]
    conv4 = [1, -2, -1, 1, -2, -2, 0, 2, 5, -1, 0, 1, -1, -1, 0, 0, 1, 2, -1, 2, -1, -1, -1, -1, 0, 2, 3, 0, -1, -1, -1, 0, 0, 1, -1, 3, -1, -3, -1, -1, -2, 2, -2, -3, 10, -1, 1, -1, -1, 2, 2, -2, -1, 4, -2, -1, 2, 1, -2, -1, 2, -2, 6, 2, 2, 2, -1, 3, -1, 3, -1, 0, 1, 2, -4, 3, 3, 2, 3, 1, 0, 2, 2, 2, 1, 1, 1, 1, -3, 0, -2, -1, -1, 2, 1, -1, 2, 1, 3, -2, 1, -2, -1, -2, 2, 2, -1, 6, 2, 2, 1, 2, -1, 1, -1, 0, 0, 1, 1, 1, -1, -2, -2, -2, -2, 7, -2, 1, 1, 1, 1, 1, 1, 1, 0, 2, 1, -2, -1, -1, -2, 1, -1, 5, -1, -1, -2, 1, 0, -2, -2, -2, 8, 3, -1, 4, 3, 2, 1, -5, 6, 0, 6, -1, 4, -1, 6, -2, -2, 5, 0, 2, 1, 2, 1, -2, 1, 1, 0, 0, -3, 1, -1, 3, 1, -2, -3, -3, 9, -3, -1, 1, -2, 1, 2, -3, -3, 9, -1, 0, -1, -1, 3, 3, 2, 3, 0]
    conv5 = [1, 0, -1, -1, 0, 0, -1, 1, 2, -1, 2, -1, -1, 0, -1, -1, 2, 3, -1, -2, 1, 1, 0, -2, -2, 2, 5, -1, -1, -1, 0, 1, 0, 0, -1, 3, -1, -1, 2, -2, 2, -1, 1, 1, 3, 2, 2, 1, 2, 1, 1, 1, -3, 0, 1, 3, 4, 4, -2, 3, -2, 1, 0, -2, 1, -1, -2, 2, 2, -2, -1, 6, 2, 1, 1, 2, 1, 1, -2, 0, 0, -2, 1, -1, 2, 1, -2, -1, -1, 5, 1, 1, -1, 1, -2, -2, -2, -2, 7, -1, -3, -1, -1, -2, 2, -2, -3, 10, 1, 1, 1, -2, 1, 1, 1, 1, 0, 2, -1, 1, -2, 2, -1, -2, -2, 6, 1, -1, -2, -3, -3, 2, 1, -3, 9, 2, 2, -1, 2, 3, -1, 3, -1, 0, -4, 2, 3, 1, 3, 2, 3, 1, 0, 4, -1, 3, 3, -5, 1, 2, 6, 0, -1, 1, -1, -1, -2, 2, 2, -1, 4, -1, 0, -1, -1, 2, 3, 3, 3, 0, 4, -1, -1, 6, -2, -2, 6, 5, 0, -1, 1, 3, -3, -3, -2, 1, -3, 9, -2, -1, 1, -1, -2, -2, 0, -2, 8]
    conv6 = [-3, -2, -3, -1, -1, -3, 1, 2, 10, -1, 0, 0, 1, -1, -1, -1, 0, 3, 1, 0, 0, -1, -1, 0, -1, -1, 3, -1, 2, 0, -1, -1, -2, -1, -2, 6, -1, -2, -3, -3, -1, 2, 1, -3, 10, -1, 2, -1, -1, 0, -2, -2, -2, 7, 1, -1, 2, 1, -2, 1, 3, 1, 0, 1, 0, -1, 1, 0, 1, 0, 1, 0, 4, 1, 1, 4, 0, -2, 3, -2, 0, 3, 4, 4, 1, 1, 3, -2, -2, 0, -1, 1, 1, -1, 1, 2, 0, 2, 0, -1, 3, -1, 2, 2, 2, 2, -1, 0, 0, 0, -1, 0, -1, 1, -1, 1, 2, 1, 4, 4, 3, 1, -2, -2, 3, 0, 2, 3, -1, -1, 2, -1, 2, 2, 0, 1, -4, 1, 1, 2, 3, 2, 3, 0, -1, 0, 1, -1, 0, -1, 1, -1, 3, 3, -4, -3, 1, 1, -3, -2, -2, 10, 1, -4, -3, 3, 1, -2, -2, -3, 10, -4, -4, 1, -1, 2, 3, -4, -2, 11, -1, -4, 1, -4, 2, -2, -4, 3, 11]
    conv7 = [-1, -2, -1, 1, -2, -1, 2, -1, 6, 2, -2, -1, -2, -1, 0, 2, 1, 4, 1, 0, 1, 2, 2, -2, 1, 2, 0, -1, -1, -2, -2, 1, -1, 2, -2, 7, -1, -2, 2, 1, 1, 1, -2, -1, 4, 3, 4, 4, 1, 2, -2, 1, -2, 0, -1, 1, 1, 1, -2, -2, -2, -1, 6, 1, -2, 1, -1, 1, -2, -2, 1, 5, 2, 1, -1, 1, 1, -2, -2, -1, 4, -2, 1, 1, -2, 2, -1, -2, 1, 5, 3, 5, 5, 4, -2, -2, -1, 1, 0, 2, 1, 2, 2, 1, 1, -3, 1, 0, 2, -1, 3, -1, 3, 2, 3, -1, 0, 2, 2, 1, 2, -3, 2, 1, 1, 0, -1, 2, -1, 3, 2, 2, 2, -1, 0, -2, -1, 0, -2, 1, 2, -1, -2, 6, 2, 1, -1, -2, -1, 2, -1, -2, 5, -2, -2, 1, 2, 1, -1, -1, -2, 6, -1, 3, -1, 3, 2, 2, -1, 2, 0, 1, -1, -1, -1, -2, 2, -1, -2, 6, -2, 2, 1, -2, -2, 2, 1, -1, 5, 0, -1, -1, -1, 3, 3, 2, 3, 0, -1, -1, -2, 1, -2, 1, -1, 2, 5, 1, -2, -2, 1, -2, -2, 1, 0, 6, -1, -1, 1, -1, -1, -1, 1, 0, 4, -1, 3, -3, -2, -1, -3, 1, 3, 7, -3, -1, -2, -2, -3, -1, -3, 1, 12]





    banzijie1 = [['0000000000000000000000000000000000000000000000000000000000000010']]
    banzijie2 = [['0000000001010000000000000000000000000000000000000000000000000000']]
                 


    imp = []


    for i in banzijie1:
        for j in banzijie2:
            PrintOuter(i, j)
            final = gubi_2(file_zuizhong)  # 第二次gurobi ，生成最小活跃S盒
            if final == GRB.Status.INFEASIBLE:
                print(i[0][::-1])
                print(j[0][::-1])



