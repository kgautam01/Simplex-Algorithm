'''
---------------------------------------------------------ASSIGNMENT-1----------------------------------------------------
                                                CS5040: LINEAR OPTIMIZATION
GROUP MEMBERS: 
KULDEEP GAUTAM      (CS20MTECH01004)
SUBODH NIGAM        (CS20MTECH01006)
VISHAL SINGH YADAV  (CS20MTECH01001)
-------------------------------------------------------------------------------------------------------------------------

QUES-2: Simplex Algorithm for the DEGENERATE case.

-------------------------------------------------------------------------------------------------------------------------
'''

import numpy as np 
def degeneracyremove(b):
    rows = len(b)
    EPSILON = 10**(-10)
    for i in range(rows):
            b[i]  = b[i] + np.random.normal() * EPSILON
    return b

def Simplex(mat_A, mat_b, mat_C):
    # A: Coeff matrix of constraints
    # b: Constant matrix of constraints (RHS)
    # C: Coeff matrix of variables in objective function
    A = np.array(mat_A).reshape(m, n)
    B = A_b = np.identity(m, dtype=int)
    b = np.array(mat_b)
    C_n = np.array(mat_C,dtype=float)
    C_b = np.array(np.zeros(m,dtype=float))

    # X: Feasible solution
    X = np.concatenate((np.zeros(n, dtype=int), b), axis=0)
    X_var = np.arange(m+n)

    while(True):
        try:
            B_inv = np.linalg.inv(B)
        except:
            print("An exception occurred of singular matrix")
            break
        
        B_inv_A = np.dot(B_inv, A)
        B_inv_b = np.around(np.dot(B_inv, b), decimals=1)
        degen = False
        if (B_inv_b == 0).any():
            degen = True
        while degen:
            B_inv_b = degeneracyremove(B_inv_b)
            if (B_inv_b == 0).any():
                degen = True
            else:
                degen = False

        B = np.dot(B_inv, B)
        A = B_inv_A
        b = B_inv_b

        # Finding feasible directions
        d_n = np.identity(n, dtype=int)
        d_b = np.negative(B_inv_A)

        # Finding improving directions
        C_n = np.add(np.dot(C_n, d_n), np.dot(C_b, d_b))
        C_b = np.zeros(m, dtype=int)

        if (C_n <= 0).all():                # Terminating Condition
            break

        else:
            e_index = np.argmax(C_n)
            pos = 0
            mini = 0
            for i in range(m):
                if d_b[i][e_index] < 0:
                    mini = -(B_inv_b[i]/d_b[i][e_index])
                    pos = i
                    break
            
            for i in range(m):
                if d_b[i][e_index] < 0:
                    ratio = -(B_inv_b[i]/d_b[i][e_index])
                    if ratio < mini:
                        mini = ratio
                        pos = i
            l_max = mini
            l_index = pos
            
            # Calculating new solution X(t+1) = X(t) + λmax*d
            d = np.concatenate((d_n[:, e_index].T, d_b[:, e_index].T), axis=0)
            if (d>=0).all():
                print("Unbounded Problem")
                break
            X = np.add(X, np.multiply(l_max, d))

            # Updating X(Feasible solution)
            X[e_index], X[n+l_index] = X[n+l_index], X[e_index]
            X_var[e_index], X_var[n+l_index] = X_var[n+l_index], X_var[e_index]

            # Updating C_n and C_b
            temp_c = np.round(C_n[e_index], decimals=0)
            C_n[e_index] = C_b[l_index]
            C_b[l_index] = temp_c

            # Updating A and B
            temp = list(B[:, l_index])
            B[:, l_index] =  list(A[:, e_index])
            A[:, e_index] = temp

    # Merging X and X_var into X_dict
    X_dict = {X_var[i]: X[i] for i in range(len(X_var))}

    # Using X_dict to find Solution Vector Z
    Z = np.zeros(n)
    for i in range(n):                  
        Z[i] = np.round(X_dict[i], decimals=2)
    if ((Z<0).any()):
        print("Infeasible Solution")
    else:                                                                                  
        print("Solution Vector: ", Z)
        print("Optimized Value: ", np.dot(np.array(mat_C), Z.T))
                    
# Taking Inputs A , b ,C
m = int(input("Enter the number of rows: "))
n = int(input("Enter the number of columns: "))
print("Enter the entries of A in a single line (separated by space): ")
mat_A = list(map(int, input().split()))
print("Enter the entries of b in a single line (separated by space): ")
mat_b = list(map(int, input().split()))
print("Enter the entries of C in a single line (separated by space): ")
mat_C = list(map(int, input().split()))
Simplex(mat_A, mat_b, mat_C)