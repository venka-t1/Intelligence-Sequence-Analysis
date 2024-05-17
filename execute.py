from main import KitBit
from input import iq_set_1,iq_set_2,sbox,kitas_sequence,accuracy_check_2,accuracy_check_1,accuracy_sbox
from math import exp

def read_path(path):
    file1 = open(path, "r")
    Content = file1.read()
    CoList = Content.split("\n")
    file1.close()
    return CoList

def write_path(path, lines):
    file1 = open(path, 'w')
    for line in lines:
        file1.write(str(line))
        file1.write('\n')
    file1.close()

def execute_model(sequences, kl, mz, path):
    results, solved = [], 0
    for seq in sequences:
        h = KitBit(seq[:-1], kl, 5000000000, 3, search_algorithm='BFS', n=2, min_zeros=mz, epsilon=exp(-18), all_solutions=False)
        x = h.solve_seq()
        results.append(x)
        if not x[0][0] or len(x[0][0])<len(seq) or x[0][0] != seq:
            continue
        else:
            solved += 1
    write_path(path, results)
    
iq1_solved_count = 0
iq1_kitas_solved_count = 0
iq1_total_count = 0
execute_model(iq_set_1, kitas_sequence, 1, '/Users/vaishnavivangireddy/Downloads/MINI_PRO_FINAL/FINAL_OUTPUT/IQ_S1Z.txt')
with open('/Users/vaishnavivangireddy/Downloads/MINI_PRO_FINAL/FINAL_OUTPUT/IQ_S1Z.txt','r') as file:
    for line, number in zip(file, accuracy_check_1):
        iq1_total_count += 1
        data = eval(line)
        if data[0][0]:
            iq1_kitas_solved_count += 1
            if data[0][0][-1] == number:
                iq1_solved_count += 1
    print("\n")
    print("IQ SEQUENCE SET - 1 \n")
    print("TOTAL : ",iq1_total_count)
    print("KITAS SOLVED COUNT : ",iq1_kitas_solved_count)
    print("PREDICTION SOLVED COUNT: ",iq1_solved_count)
    accuracy = (iq1_solved_count/ iq1_total_count) * 100
    print("PERCENTAGE OF KITA SEQUENCE GENERATION : ", (iq1_kitas_solved_count/iq1_total_count)*100)
    print("ACCURACY OF NUMBER PREDICTION :", accuracy)
    print("---------------------")
    print("\n")

iq2_total_count,iq2_kitas_solved_count,iq2_solved_count=0,0,0
execute_model(iq_set_2, kitas_sequence, 1, '/Users/vaishnavivangireddy/Downloads/MINI_PRO_FINAL/FINAL_OUTPUT/LI_S1Z.txt')
with open('/Users/vaishnavivangireddy/Downloads/MINI_PRO_FINAL/FINAL_OUTPUT/LI_S1Z.txt','r') as file:
    for line, number in zip(file, accuracy_check_2):
        iq2_total_count += 1
        data = eval(line)
        if data[0][0]:
            iq2_kitas_solved_count += 1
            if data[0][0][-1] == number:
                iq2_solved_count += 1
    print("\n")
    print("IQ SEQUENCE SET - 2 \n")
    print("TOTAL : ",iq2_total_count)
    print("KITAS SOLVED COUNT : ",iq2_kitas_solved_count)
    print("PREDICTION SOLVED COUNT : ",iq2_solved_count)
    accuracy = (iq2_solved_count/ iq2_total_count) * 100
    print("PERCENTAGE OF KITA SEQUENCE GENERATION : ", (iq2_kitas_solved_count/iq2_total_count)*100)
    print("ACCURACY OF NUMBER PREDICTION :", accuracy)
    print("-----------------------")
    print("\n")
total_count,kitas_solved_count,solved_count=0,0,0
execute_model(sbox, kitas_sequence, 1, '/Users/vaishnavivangireddy/Downloads/MINI_PRO_FINAL/FINAL_OUTPUT/SBOX.txt')
with open('/Users/vaishnavivangireddy/Downloads/MINI_PRO_FINAL/FINAL_OUTPUT/SBOX.txt','r') as file:
    for line, number in zip(file, accuracy_sbox):
        total_count += 1
        data = eval(line)
        if data[0][0]:
            kitas_solved_count += 1
            if data[0][0][-1] == number:
                solved_count += 1
    print("\n")
    print("DES S-BOX PREDICTION \n")
    print("TOTAL : ",total_count)
    print("KITAS SOLVED COUNT : ",kitas_solved_count)
    print("PREDICTION SOLVED COUNT : ",solved_count)
    accuracy = (solved_count/ total_count) * 100
    print("PERCENTAGE OF KITA SEQUENCE GENERATION : ", (kitas_solved_count/total_count)*100)
    print("ACCURACY OF NUMBER PREDICTION :", accuracy)
    print("-----------------------")
    print("\n")

def execute_OEIS(sols, kl, sols_cad):
    solved, not_solved, q = [], [], 0
    for i in range(len(sols)):
        if i % 100 == 0:
            print("SEQUENCE NO : "+str(i)+"   "+"SOLVED: "+str(len(solved))+"   "+"NOT SOL: "+str(len(not_solved)))
        seq, solc = sols[i], sols_cad[i]
        kitasi = kl[i][1:-1].split('&')
        h = KitBit(seq[:-1], kitasi, 5000000000, kitasi, search_algorithm='branch', n=1, min_zeros=1, epsilon=exp(-18), all_solutions=False)
        x = h.solve_seq()
        if not x[0][0] or len(x[0][0])<len(seq) or x[0][0] != seq:
            h = KitBit(seq[:-1], kitasi, 5000000000, kitasi, search_algorithm='branch', n=1, min_zeros=2, epsilon=exp(-18), all_solutions=False)
            y = h.solve_seq()
            if not y[0][0] or y[0][0] != seq or len(y[0][0])<len(seq):
                not_solved.append(solc)
            else:
                q += 1
                solved.append(y)
        else:
            solved.append(x)
    print("PERCENTAGE OF KITAS SEQUENCE PREDICTION : ",((len(solved)/(len(solved)+len(not_solved)))*100))
    return solved, not_solved

sols, sols_cad = [], []
OEIS_data = read_path('/Users/vaishnavivangireddy/Downloads/MINI_PRO_FINAL/FINAL_INPUT/OEIS_SERIES_SOLVED.txt')
OEIS_kitas = read_path('/Users/vaishnavivangireddy/Downloads/MINI_PRO_FINAL/FINAL_INPUT/OEIS_KITAS.txt')
print("OEIS SEQUENCES : ",len(OEIS_data))
print("---------------------------------------------------")
for j in range(len(OEIS_data)):
    seq1 = list(map(lambda u: int(u), OEIS_data[j][9:-2].split(',')))
    pos_sol = seq1[:]
    if len(pos_sol) < 4 or len(seq1[:-1]) < 2:
        continue
    sols.append(pos_sol)
    sols_cad.append(OEIS_data[j][:-1])

sol_def, not_sol = execute_OEIS(sols, OEIS_kitas, sols_cad)
write_path('/Users/vaishnavivangireddy/Downloads/MINI_PRO_FINAL/FINAL_OUTPUT/OEIS_results.txt', sol_def)