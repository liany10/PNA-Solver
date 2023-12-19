# coding:utf-8
import re


class sequence(object):

    def __init__(self, st):
        self.seq = st
        self.seqType = "Unknown"
        
##    def __init__(self, st, seqType):
##        self.seq = st
##        self.seqType = seqType

    def __str__(self):
        return self.seq
    def isPNA(self):        
        for monomer in ('L', 'l','Q','q','O','o','R','r','S','s','E','e'):
            if self.seq.find(monomer) != -1:
                return true
        return false
    
    def isDNA(self):
        if self.isPNA():
            return false
        return (self.seq.find('t') != -1) or (self.seq.find('T') != -1)
    def isRNA(self):
        return not self.isDNA()
        

    def Reverse(self):
        return self.seq[::-1]

    def DnaToRna(self):
        return self.rep('t', 'u')

    def RnaToDna(self):
        return self.rep('u', 't')

    def Complementary(self):
        st = ""
        for i in self.seq:
            st = st + self.basePair(i)
        return st

    def basePair(self, b):
        if b == 'a' and self.isDNA():
            return 't'
        if b == 'a' and self.isRNA():
            return 'u'        
        if b == 't':
            return 'a'
        if b == 'u':
            return 'a'        
        if b == 'c':
            return 'g'
        if b == 'g':
            return 'c'
        if b == 'A' and self.isDNA():
            return 'T'
        if b == 'A' and self.isRNA():
            return 'U'        
        if b == 'T':
            return 'A'
        if b == 'U':
            return 'A'        
        if b == 'C':
            return 'G'
        if b == 'G':
            return 'C'        
        return 'X'
    

    def rep(self, a, b):    
        return self.seq.replace(a, b).replace(a.upper(), b.upper())

    def extinctionCoefficientCalc(self):  ##s2U用R代替，S指Ssuc  数值参照https://www.novopro.cn/tools/oligo-calculation.html 默认字典是针对PNA的
        coefficient260nmDict = {'A':15.4, 'C':7.3, 'G':11.7, 'T':8.8, 'O':9.4, 'R':10.2,'2':10.2, 'S':16.4, 'L':7.3, 'Q':7.3, 'E':6.0}
        coefficient260nmDictRNA = {'A':15.4, 'C':7.2, 'G':11.5, 'U':9.9}
        coefficient260nmDictDNA = {'A':15.4, 'C':7.4, 'G':11.5, 'T':8.7}
        tmpseq = self.seq.upper()   ##统一变大写
        st = 0
        for i in tmpseq:
            st = st + coefficient260nmDict[i]
        return ((float)(st))
    
    
    def MWCalc(self):  ##s2U用R代替，S指Ssuc
        coefficientMWDict = {'A':275.27, 'C':251.25, 'G':291.27, 'T':266.26, 'O':9.4, 'R':10.2,'2':10.2, 'S':16.4, 'L':267.31, 'Q':350.38, 'E':6.0}
        tmpseq = self.seq.upper()   ##统一变大写
        st = 0
        for i in tmpseq:
            st = st + coefficientMWDict[i]
        return ((float)(st +145.2))    

    def ConcentrationCalc(self, ValueA260):
        return ((float)(ValueA260/self.extinctionCoefficientCalc()))
