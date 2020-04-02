import re
from Bio.Seq import translate,Seq

li = open('yarrowia1.fa').readlines()

start_line = li[0]
length = len(li)
num = []
for i in range(1, length - 1):
    if start_line.split('\t')[0] != li[i].split('\t')[0]:
        num.append(i)
        start_line = li[i]

num.insert(0, 0)
num.append(length - 1)
num.append(length)

total = []
for i in range(len(num) - 1):
    s = ''
    l = ''
    ls = [0]
    n = 0
    for j in (li[num[i]:num[i + 1]]):
        s += j.split('\t')[5].rstrip('\n')
        ls.append(n +len(j.split('\t')[5].rstrip('\n')))
        l = j.split('\t')[0]
        n += len(j.split('\t')[5].rstrip('\n'))
    total.append([l, s, str(translate(Seq(s))), ls])
    #print('\n')

#print(total) 
t = []
for i in total:
    index = []
    if 'Q' in i[2]:
        index = [m.start()*3 for m in re.finditer("Q", i[2])]
#        print(index)
#        index.append(i[2].index('Q'))
#    else:
#        continue
        t.append([i[0], i[1], i[2], i[3],index])
#print(t)

to1 = []
for i in t:
    #print(i[4][0])
#    print(i)
    d1 = {}
#    cr = []
    for j in i[4]:
        cr = []
        for n in range(len(i[3]) - 1):            
            if i[3][n] <= j -10 and j + 13 <= i[3][n + 1]:
                seq = (i[1][j -10:j + 13].upper())
                if seq[0 : 2] == 'TT' and seq[10:16].count('C') == 1 and seq[9] != 'G':
                    cr.append(seq)
        d1.update({j:cr})
    to1.append([i[0], i[1], i[2], i[3], i[4], d1])                
                
#print(to1)
#for i in to1:
#    print(i)


to2 = []
for i in t:
    d1 = {}
#    cr = []
    for j in i[4]:
        cr = []
        for n in range(len(i[3]) - 1):
            if i[3][n] <= j -11 and j + 12 <= i[3][n + 1]:
                seq = (i[1][j -11:j + 12].upper())
                if seq[0 : 2] == 'TT' and seq[10:16].count('C') == 1 and seq[10] != 'G':
                    cr.append(seq)
#                print(seq)
        d1.update({j:cr})
    to2.append([i[0], i[1], i[2], i[3], i[4], d1])
#print(to2) 
#for i in to2:
#    print(i)

to3 = []
for i in t:
    d1 = {}
    for j in i[4]:
        cr = []
        for n in range(len(i[3]) - 1):
            if i[3][n] <= j -12 and j + 11 <= i[3][n + 1]:
                seq = (i[1][j -12:j + 11].upper())
                if seq[0 : 2] == 'TT' and seq[10:16].count('C') == 1 and seq[11] != 'G':
                    cr.append(seq)
        d1.update({j:cr})
    to3.append([i[0], i[1], i[2], i[3], i[4], d1])

#for i in to3:
#    print(i)
to4 = []
for i in t:
    d1 = {}
    for j in i[4]:
        cr = []
        for n in range(len(i[3]) - 1):
            if i[3][n] <= j -13 and j + 10 <= i[3][n + 1]:
                seq = (i[1][j -13:j + 10].upper())
                if seq[0 : 2] == 'TT' and seq[10:16].count('C') == 1 and seq[12] != 'G':
                    cr.append(seq)
        d1.update({j:cr})
    to4.append([i[0], i[1], i[2], i[3], i[4], d1])


def func(l): 
    for i in l:
        for j in (i[5]):
        #print(j)    
            if i[5][j] != []:
                print(str(i[0]) + '\t' + str(i[1]) + '\t' + str(i[2]) + '\t' + str(i[3]) + '\t' + str(i[4]) + '\t' +
                       str(j) + '\t' + i[5][j][0])
#    print(i)
func(to1)
func(to2)
func(to3)
func(to4)

to5 = []
for i in t:
    d1 = {}
    for j in i[4]:
        cr = []
        for n in range(len(i[3]) - 1):
            if i[3][n] <= j -14 and j + 11 <= i[3][n + 1]:
                seq = (i[1][j -14:j + 9].upper())
                if seq[0 : 2] == 'TT' and seq[10:16].count('C') == 1 and seq[13] != 'G':
                    cr.append(seq)
        d1.update({j:cr})
    to5.append([i[0], i[1], i[2], i[3], i[4], d1])

func(to5)
to6 = []
for i in t:
    d1 = {}
    for j in i[4]:
        cr = []
        for n in range(len(i[3]) - 1):
            if i[3][n] <= j -15 and j + 8 <= i[3][n + 1]:
                seq = (i[1][j -15:j + 8].upper())
                if seq[0 : 2] == 'TT' and seq[10:16].count('C') == 1 and seq[14] != 'G':
                    cr.append(seq)
        d1.update({j:cr})
    to6.append([i[0], i[1], i[2], i[3], i[4], d1])

func(to6)
