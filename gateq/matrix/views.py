from django.http import HttpRequest, HttpResponse, QueryDict
from django.shortcuts import render
import numpy as np
import math

# Create your views here.
ket0=np.array([1,0]);
ket1 =np.array([1,0]);

ket0=np.resize(ket0,(2,1))
ket1=np.resize(ket1,(2,1))

ket00=np.kron(ket0,ket0)
ket000=np.kron(ket00,ket0)
SX = np.array([[0,1],[1,0]]) #Gate1 2x2
I=np.array([[1,0],[0,1]])    #Gate2 2x2

H=np.array([[1,1],[1,-1]]) 

Cxx = np.kron(np.dot(ket0,np.matrix.getH(ket0)), np.kron(I,I)) + np.kron(np.kron(np.dot(ket1,np.matrix.getH(ket1)),SX),SX) #Gate5. 8x8
Cnot= np.kron(np.dot(ket0,np.matrix.getH(ket0)),np.kron(I,I)) + np.kron(np.kron(np.dot(ket1,np.matrix.getH(ket1)),I),SX) #1,3.  #Gate6 8x8


#Rotation Gate PG1 Gate7 2x2
def R(thi):
    thi=float(thi)
    a= 2 * thi;
    R=np.array([[math.cos(a), 1j*math.sin(a)], [1j*math.sin(a), math.cos(a)]]);
    return R

#C-x,y gate PG2 Gate8 8x8
def cxy(pi):
    pi=float(pi)
    th2=(math.pi)/pi
    Rxy=np.array([[math.cos(th2/2), 0, 0, 1j*math.sin(th2/2)],[0, math.cos(th2/2), 1j*math.sin(th2/2), 0], [0, 1j*math.sin(th2/2), math.cos(th2/2), 0], [1j*math.sin(th2/2), 0, 0, math.cos(th2/2)]])
    Cxy= np.kron(np.dot(ket0,np.matrix.getH(ket0)),np.kron(I,I)) + np.kron(np.dot(ket1,np.matrix.getH(ket1)),Rxy)
    return Cxy

def index(request):
    parameters=''
    if request.method == 'POST':
        a=request.POST.get('21')
        data=request.POST.dict()
        print(data)
        inputs= getlist(data)
        print('inputs',inputs)
        product,states,values=Quantum_Operator(inputs)
        print(product,states,product.shape)
        parameters=zip(product,states,values)
    return render(request,'index (1).html',{'parameters':parameters})

def getlist(data):
    l1=list(data.values())
    input_list=l1
    inputs=input_list[:-1]
    return inputs

def retrieve_matrix(gate):
    if gate=='H':
        return H
    if gate=='I':
        return I
    if gate=='SX':
        return SX
    
    if 'Cxx' in gate:
        
        return Cxx
    
    if 'R' in gate:
        stop=gate.find(')')
        thi=gate[2:stop]
        return R(thi)
    
    if 'Cxy' in gate:
        stop=gate.find(')')
        pi=gate[4:stop]
        return cxy(pi)
    
    
    
    if gate=='Cnot':
        return Cnot
    
    
def Kr_prod(kr_input):
    
    x=np.kron(retrieve_matrix(kr_input[0]),retrieve_matrix(kr_input[1]))
    y=np.kron(x,retrieve_matrix(kr_input[2]))
    
    return y

# Matrix Multiplication
def Matrix_Mult(x,y):
    mult=np.dot(x,y)
    return mult


def Quantum_Operator(inputs):
    
    
        input_list=inputs

        multK=[]
        K_list=[]
        num=0
        
        for i in range(len(inputs)):
            K_list.append(inputs[i])
            num=num+1
            if num%3==0:
                multK.append(K_list)
                K_list=[]
    
       
        Final_Matrices=[]

        n_col=len(multK)

        for i in range(n_col):
            
            if 'Cnot' in multK[i]:
                final=retrieve_matrix('Cnot')
            elif 'Cxx' in multK[i]:
                final=retrieve_matrix('Cxx')
            
            elif 'Cxy' in multK[i][0]:
                final=retrieve_matrix(multK[i][0])
            else:
                final=Kr_prod(multK[i])
    
            Final_Matrices.append(final)
        
        Final_Matrices=list(reversed(Final_Matrices))
        
        from functools import reduce

        product=reduce(Matrix_Mult,Final_Matrices)
        
        product=np.dot(product,ket000)
        
        ####Changes####
        
        norm=np. linalg. norm(product)
        final=product/norm
        index=0
        
        states=['000','001','010','011','100','101','110','111']
        values=list(np.diag(final.dot(np.conjugate(final.T))))
        list1=[0,0,0,0,0,0,0,0]
        for n,ele in enumerate(values):
            print(f"|{states[index]}>:P[{str(index)}]={abs(ele)}")
            a=str('|'+states[index]+'>:P['+str(index)+']='+str(abs(ele)))
            list1[n]=a
            index=index+1
        # import matplotlib.pyplot as plt
        # states=['000','001','010','011','100','101','110','111']
        # plt.xlabel('States')
        # plt.ylabel('Probability Score')
        # plt.bar(states,values)
        
        return product,states,values
