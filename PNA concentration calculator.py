# coding:utf-8
import re
import sequence
##import RNA
from subprocess import Popen, PIPE

class XNADealer():

    def __init__(self):
        self.RNAs = {}  ## rna名字都带">"!
        self.PNAs = {}
        self.RNAStructures = {}
        
    def SetRNAsAndPNAs(self, rnatxtpath, pnatxtpath):
        self.RNAs = self.ReadFASTAIntoDict(rnatxtpath)
        self.PNAs = self.ReadFASTAIntoDict(pnatxtpath) ## from N to C, use α for asPNA and β for dbPNA
        for pna in self.PNAs:
            self.PNAs[pna] = self.PNAs[pna].replace("伪","α").replace("尾","β")
##        FindExactPairs(PNAs, premiRNAs, "outputs/results2.txt")        

    def SetRNAStructures(self, txtpath):
        self.RNAStructures = self.ReadFASTAIntoDict(txtpath)
        
    def ReadPNAandCalculate(self):
        while True:
            seq1 = sequence.sequence(input("Input sequence:"))
            if (seq1.seq == "") or (seq1.seq == "N"):
                break
            print("Extinction coefficient:" + str(seq1.extinctionCoefficientCalc()))
            print("MW(only ATCGLQ allowed):" + str(seq1.MWCalc()))

    def ReadFASTAIntoDict(self, filename):
        XNAs = {}
        name = ""
        sequence = ""
        f = open(filename,'r')
        while True:
            line = f.readline()
            if not line:
                break
            if (line[0] == ">"):
                name = line.strip()
                sequence = ""
            else:       ##将序列第二行连接上
                sequence += line.strip()
            XNAs[name] = sequence
        f.close()
        return XNAs

    def asPNAFowardComplementary(self, asPNAseq):
        asPNAFowardPair = {'T':'A', 'C':'G', 'A':'U', 'G':'C', '2':'A'} ##加了S2U，不知是否合适
        st = ""
        for base in asPNAseq:
            if base in asPNAFowardPair:
                st = st + asPNAFowardPair[base]
            else:
                print("Wrong base in asPNA part!")
            ##处理st
        return st

    def asPNAReverseComplementary(self, asPNAseq):
        return self.asPNAFowardComplementary(asPNAseq)[::-1]

    def dbPNAFowardComplementary(self, dbPNAseq):
        dbPNAFowardPair = {'L':'G', 'Q':'C', 'T':'A', '2':'A', 'R':'A'} ## T,2,R一样
        st = ""
        for base in dbPNAseq:
            if base in dbPNAFowardPair:
                st = st + dbPNAFowardPair[base]
            else:
                print("Wrong base in dbPNA part!")
        return st
    def dbPNAFowardComplementaryOtherSide(self, dbPNAseq):
        dbPNAFowardPairOtherSide = {'L':'C', 'Q':'G', 'T':'U', '2':'U', 'R':'U'} ## T,2,R一样
        st = ""
        for base in dbPNAseq:
            if base in dbPNAFowardPairOtherSide:
                st = st + dbPNAFowardPairOtherSide[base]
            else:
                print("Wrong base in dbPNA part!")
        return st

    def dbPNAReverseComplementary(self, dbPNAseq):
        return self.dbPNAFowardComplementary(dbPNAseq)[::-1]
    def dbPNAReverseComplementaryOtherSide(self, dbPNAseq):
        return self.dbPNAFowardComplementaryOtherSide(dbPNAseq)[::-1]

    def pairSequences(self, seqclass): ##return a tuple PNA
        seq = seqclass.seq
        a = seq.find("α")
        b = seq.find("β")
        if a < b: ## asPNA在dbPNA前(N到C) e.g. αCAACAGβTQTTQT
            asPNA = seq[seq.find("α") + 1:seq.find("β")]
            dbPNA = seq[seq.find("β") + 1:]
            return((self.dbPNAFowardComplementary(dbPNA), self.dbPNAReverseComplementaryOtherSide(dbPNA) + self.asPNAReverseComplementary(asPNA)))
        else:
            dbPNA = seq[seq.find("β") + 1:seq.find("α")]
            asPNA = seq[seq.find("α") + 1:]
            return((self.dbPNAFowardComplementary(dbPNA),self.asPNAReverseComplementary(asPNA) + self.dbPNAReverseComplementaryOtherSide(dbPNA))) 

    def PrintAndWrite(self, st, file):
        print(st)
        file.write(st + "\n")

    def FindExactPairsInSpecies(self, filename, species):       ##可改为正则条件; 目前未考虑preRNA结构
        fw = open(filename, "w")
        
        for pna in self.PNAs:
            self.PrintAndWrite(pna + ":" +self.PNAs[pna], fw)
            
            PNAab = self.pairSequences(sequence.sequence(self.PNAs[pna]))
            self.PrintAndWrite("PNA pairing: " + PNAab[0] + "  " + PNAab[1], fw)
            for premiRNA in self.RNAs:
                if premiRNA.find(species) == -1:
                    continue
                premiRsequence = self.RNAs[premiRNA]            
                if premiRsequence.find(PNAab[0], 0, len(premiRsequence)) != -1 and premiRsequence.find(PNAab[1], 0, len(premiRsequence)) != -1:
                    self.PrintAndWrite(premiRNA + ":\n" + premiRsequence, fw)
                    
        fw.close()

    def CreateAndWriteRNAStructure(self, filename):
        fw = open(filename, "w")
        
        for rna in self.RNAs:            
            if rna.find("hsa") == -1:
                continue
            fw.write(rna + "\n")
            p = Popen('RNAfold/RNAfold.exe', stdin=PIPE, stdout=PIPE)
            answer = p.communicate(self.RNAs[rna].encode())          
            ans = answer[0].decode('utf-8')
            fw.write(ans[ans.find('\n')+1: ans.find(' ')] + "\n")
        fw.close()
        
    def FindRNAsFordbPNA(self, length, bulgeallowed):
        ret = []
        for rna in self.RNAStructures:
            ans53 = self.LongestdbPNATarget(self.RNAs[rna], self.RNAStructures[rna], bulgeallowed, "(","AG")
            ans35 = self.LongestdbPNATarget(self.RNAs[rna], self.RNAStructures[rna], bulgeallowed, ")","AG")
            if ans53[0] >= ans35[0]: ##将较长序列作为最终答案
                ans = ans53
            else:
                ans = ans35
            if ans[0]>= length:
                ret.append(rna)
                ret.append(self.RNAs[rna])
                ret.append(self.RNAStructures[rna])
                ret.append(ans)
                ret.append(self.RNAs[rna][ans[1]:ans[2]+1])
        return ret

    def LongestdbPNATarget(self, seq, structure, bulgeallowed, allowedBase, allowedStructure): ## 返回三元tuple(序列长度，起点，终点),只允许"("; 最长长度
        cur_start = -1
        
        for i in range(len(seq)):
            if self.SuitsCondition(seq, structure, i, allowedBase,allowedStructure):
                cur_start = i
                cur_end = i
                cur_len = cur_end - cur_start + 1
                ans_start = cur_start
                ans_end = cur_end
                ans_len = cur_len
                break
        if cur_start == -1:
            return(0, -1, -1)
        i = i + 1
        while i < len(seq):            
            while self.SuitsCondition(seq, structure, i, allowedBase,allowedStructure):    ## rna[i]符合条件, i++
                cur_end = i
                cur_len = cur_end - cur_start + 1
                i = i + 1
            if ans_len < cur_len:           ## 更新ans
                ans_start = cur_start
                ans_end = cur_end
                ans_len = cur_len
            while (i < len(seq)) and (not self.SuitsCondition(seq, structure, i, allowedBase,allowedStructure)): ## i++, 越过不符条件区域
                i = i + 1
            if i >= len(seq):    ##整条序列5'>3'方向检查完毕
                if bulgeallowed:
                    return (ans_len, ans_start, ans_end)
                else:
                    a = self.BulgePosFor(seq, structure, ans_len, ans_start, ans_end)
                    if a == -1:
                        return (ans_len, ans_start, ans_end)
                    else:
                        return(a, ans_start, ans_start - 1 + a *((ans_end-ans_start)//abs(ans_end-ans_start))) ##返回第一个bulge前的长度. 此处导致e.g.第一个碱基就是bulge时无法返回更好的结果
            cur_start = i
        

    def BulgePosFor(self, seq, structure, le, st, en): ##检查rna的st到en的序列配对区域是否含'.',即bulge结构, 返回第一个bulge前的长度，无则返回-1
        seqBracket = structure[st]## 待检查原序列方向是"("还是")"
        if (seqBracket == '('):
            direction = 1
            i = en + direction
        else:
            direction = -1
            i = st + direction
        
        midpairtmp = 0
        while ((i >= 0) and (i < len(structure))): ##找到配对序列起始坐标
            if (structure[i] == '.'):
                i = i + direction
            elif (structure[i] == seqBracket):  ##碰到同向括号，将之加入待配对数
                midpairtmp = midpairtmp + 1
                i = i + direction
            elif (midpairtmp > 0):  ##异向括号，但还剩括号没配对
                midpairtmp = midpairtmp -1
                i = i + direction
            else:  ##异向括号，且中间序列全配对了
                break
            
        for k in range(0, le):
            if structure[i + k * direction] == '.':
                return k
        return -1
        
    def SuitsCondition(self, seq, structure, index, structurecondition, rnaseqcondition):
        if (index < 0) or (index >= len(seq)):
            return False
        return (structure[index] in structurecondition) and (seq[index] in rnaseqcondition)
        
if __name__ == "__main__":
    dealer = XNADealer()
    dealer.SetRNAsAndPNAs("inputs/all miRNA hairpins.fa.txt", "inputs/PNAs.txt")
    dealer.SetRNAStructures("inputs/hsa miRNA hairpin structures.txt")
    
##    dealer.ReadPNAandCalculate()
    
##    dealer.FindExactPairsInSpecies("outputs/results3.txt", "hsa")
##    dealer.CreateAndWriteRNAStructure("outputs/all miRNA hairpin structures.txt")
    
    
##    ans = dealer.FindRNAsFordbPNA(12, False)
##    for a in ans:
##        print(a)
    a = sequence.sequence("----·-··---··--- -----·----- -··-  ")
    print(a.Reverse())


##    if len(ans)>0:
##        for i in range(0,5):
##            print(ans[i])
            
            

    
