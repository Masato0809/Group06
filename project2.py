# import libraries
import chemparse, fractions, math

def balancing(txt):
    # divide into left and right sides
    reactants, products = [elem.replace(" ","") for elem in txt.split("=")]

    # parse the chemical formulas using chemparse
    reactant_compounds = [chemparse.parse_formula(comp) for comp in reactants.split("+")]
    product_compounds = [chemparse.parse_formula(comp) for comp in products.split("+")]

    # extract unique elements
    elements = sorted(set().union(*[comp.keys() for comp in reactant_compounds + product_compounds]))

    elements_1 = sorted(set().union(*[comp.keys() for comp in reactant_compounds]))
    elements_2 = sorted(set().union(*[comp.keys() for comp in product_compounds]))

    # check if the input equation does exist
    if elements_1 == elements_2:
        pass
    else:
        print("The equation does not exist!")
        exit()

    # create matrix
    aug_matrix = []
    for elem in elements:
        row = []
        for comp in reactant_compounds:
            row.append(comp.get(elem, 0))
        for i, comp in enumerate(product_compounds):
            row.append(-1*comp.get(elem, 0))
        aug_matrix.append(row)

    # check if the rows in the matrix are not identical
    remove_row = []
    for i in range(len(aug_matrix)):
        for j in range(i+1, len(aug_matrix)):
            if aug_matrix[i] == aug_matrix[j]:
                remove_row.append(j)
            for k in range(2, 6):
                if aug_matrix[i] == k*aug_matrix[j]:
                    remove_row.append(j)
    for k in remove_row:
        aug_matrix.pop(k)

    # determine the length of rows and columns of the augmented matrix
    n_rows, n_cols = len(aug_matrix), len(aug_matrix[0])

    # detect if the augmented matrix is created
    if n_rows + 2 <= n_cols:
        print("This equation cannot be balanced by this program!")
        exit()

    # gauss-jordan elimination
    aug_matrix = aug_matrix[:-(n_rows-n_cols+1)] if n_rows > n_cols - 1 else aug_matrix
    for i in range(n_cols-1):
        # find pivot row
        # a lambda function defines a single argument "j" in the range()
        pivot_row = max(range(i, n_cols-1), key = lambda j: abs(aug_matrix[j][i]))

        # swap the pivot row with the current row
        aug_matrix[i], aug_matrix[pivot_row] = aug_matrix[pivot_row], aug_matrix[i]

        # check for zero pivot element
        pivot = aug_matrix[i][i]
        if pivot == 0:
            print("This equation cannot be balanced by this program!")
            exit()

        # scale pivot row
        aug_matrix[i] = [elem * 1/pivot for elem in aug_matrix[i]]

        # eliminate other rows
        for j in range(n_cols-1):
            if i != j:
                scale = -aug_matrix[j][i]
                aug_matrix[j] = [elem1 + scale * elem2 for elem1, elem2 in zip(aug_matrix[j], aug_matrix[i])]

    # compute coefficients for each compound
    coeffs = []
    for i in range (len(aug_matrix)):
        if i <= len(reactants)-1:
            coeffs.append(aug_matrix[i][-1]*-1)
    coeffs.append(1)

    # covert coeffs to fractions
    coeffs_1 = []
    for i in range(n_cols):
        coeffs_1.append(fractions.Fraction(coeffs[i]).limit_denominator(10**6))

    # compute LCM of coefficients
    lcm_coeff = coeffs_1[0].denominator
    for i in range(1, n_cols):
        if i == 1:
            a, b = lcm_coeff, coeffs_1[i].denominator
            a = math.lcm(a, b)
        else:
            b = coeffs_1[i].denominator
            a = math.lcm(a, b)
    lcm = a

    # multiply coefficients by LCM
    coeffs = [int(coeff*lcm) for coeff in coeffs_1]

    # generate a balanced equation string
    len_i = len(reactants.split("+"))
    balanced_equation = ''
    for i, compound in enumerate(reactants.split('+')):
        balanced_equation += f'{coeffs[i] if coeffs[i]!=1 else ""}{compound}'
        if i != len_i - 1:
            balanced_equation += ' + '
    balanced_equation += ' = '
    for i, compound in enumerate(products.split('+')):
        balanced_equation += f'{coeffs[i+len_i] if coeffs[i+len_i]!=1 else ""}{compound}'
        if i != len(products.split('+')) - 1:
            balanced_equation += ' + '

    # print the balanced equation
    print(f"Balanced equation: {balanced_equation}")

    return balanced_equation
    
def is_combustion(equation, balanced_equation):
    # split equation into reactants and products
    reactants, products = equation.split('=')
    
    # conditions to detect combustion reaction 
    O2_reactant = "O2" in reactants
    C_reactant = "C" in reactants
    H2O_product = "H2O" in products
   
    CO_product = "CO" in products
    C_product = "C" in products 

    # conditions to check if its complete or incomplete
    if balanced_equation == "2H2 + O2 = 2H2O":
        print("This is a combustion reaction equation")
        print("The reaction is complete")
    elif O2_reactant and C_reactant and H2O_product:
        print ("This is a combustion reaction equation")
        if CO_product or C_product:
            print("The reaction is complete")
        else:
            print("The reaction is incomplete")
    else:
        print("This is not a combustion reaction equation")

def find_pH_1():
    print("Units must be consistant!")
    Ma = float(input("Enter molar concentration of a strong acid: "))
    Va = float(input("Enter volume of a strong acid: "))
    Mb = float(input("Enter molar concentration of a strong base: "))
    Vb = float(input("Enter volume of a strong base: "))
    pH = 7-math.asinh((Ma*Va-Mb*Vb)/(2*math.sqrt(10**(-14))*(Vb+Va)))/math.log(10)
    print(f"The pH is {pH}")

def find_pH_2():
    print("Units must be consistant!")
    Mb = float(input("Enter molar concentration of a strong base: "))
    Vb = float(input("Enter volume of a strong base: "))
    Ma = float(input("Enter molar concentration of a strong acid: "))
    Va = float(input("Enter volume of a strong acid: "))
    pH = 7+math.asinh((Mb*Vb-Ma*Va)/(2*math.sqrt(10**(-14))*(Va+Vb)))/math.log(10)
    print(f"The pH is {pH}")

def find_pH_3():
    print("Units must be consistant!")
    Ma = float(input("Enter molar concentration of a weak monoprotic acid: "))
    Va = float(input("Enter volume of a weak monoprotic acid: "))
    Mb = float(input("Enter molar concentration of a strong base: "))
    Vb = float(input("Enter volume of a strong base: "))
    Veq = float((Ma*Va)/Mb)
    Ka = float(input("Enter the acid dissociation constant: "))
    β = (((((10**(-28))*((Va+Veq)**2))/(2*(Mb**2)*((Vb+Ka*(Va+Veq)/Mb)**2)))*(1+math.sqrt(1+4*((Ka*Mb*(Veq+(Ka*(Va+Veq)/Mb)))/((10**(-14))*(Va+Veq))))))+(Ka*(10**(-14))*(Va+Veq))/(Mb*(Veq+(Ka*(Va+Veq)/Mb))))**(-1/2)
    pH = (math.log10(β))-(math.asinh((((Veq-Vb)*Mb/(Va+Vb))+((10**(-14))/Ka))/2*β*(10**(-14)))/math.log(10))
    print(f"The pH is {pH}")

def find_pH_4():
    print("Units must be consistant!")
    Mb = float(input("Enter molar concentration of a weak monoprotic base: "))
    Vb = float(input("Enter volume of a weak monoprotic base: "))
    Ma = float(input("Enter molar concentration of a strong acid: "))
    Va = float(input("Enter volume of a strong acid: "))
    Veq = float((Mb*Vb)/Ma)
    Kb = float(input("Enter the base dissociation constant: "))
    β = ((((10**(-28))*((Vb+Veq)**2)/2*(Ma**2)((Va+Kb*(Vb+Veq)/Ma)**2))*(1+math.sqrt(1+4*(Kb*Ma(Veq+Kb(Vb+Veq)/Ma))/(10**(-14))*(Vb+Veq))))+(Kb*(10**(-14))(Vb+Veq))/Ma(Veq+Kb(Vb+Veq)/Ma))**(-1/2)
    pH = -math.log10(β)-math.asinh((((Veq-Va)*Ma/(Vb+Va))+((10**(-14))/Kb))/2*β*(10**(-14)))/math.log(10)
    print(f"The pH is {pH}")

# ask whether a balancing or pH
print("1 - balancing a equation and check if it is combustion \n2 - find out a pH")

choice = input("Which do you choose? (1/2): ")
while choice not in ["1", "2"]:
    choice = input("Please enter either 1 or 2: ")

# balance an input equation and detect combustion reaction
if choice == "1":
    # ask the user to input a chemical equation with an equal sign
    while True:
        txt = input("Enter a chemical equation without any coefficients: ")
        if "=" not in txt:
            print("Error: the equation must contain an equal sign!")
        else:
            break

    balanced_equation = balancing(txt)
    is_combustion(txt, balanced_equation)

# calculate pH
elif choice == "2":
    print("1: Titration of a strong acid with a strong base\n2: Titration of a strong base with a strong acid")
    print("3: Titration of a weak monoprotic acid with a strong base\n4: Titration of a weak monoprotic base with a strong acid")
    choice = int(input("Which type of titration would you like to deal with? (1/2/3/4): "))
    if choice == 1:
        find_pH_1()
    elif choice == 2:
        find_pH_2()
    elif choice == 3:
        find_pH_3()
    elif choice == 4:
        find_pH_4()
    else:
        print("Invalid input")

else:
    print("Invalid input")


"""
balanced_equation = balancing("Na3PO4 + MgCl2 = NaCl + Mg3(PO4)2")
is_combustion("Na3PO4 + MgCl2 = NaCl + Mg3(PO4)2", balanced_equation)
"""
# C2O2H2+HCHO+NH3=H2O+C3N2H4
# C3H8O+C13H19O8I=C3H6O+C9H9O4I+C2H6O2
# Na3PO4+MgCl2=NaCl+Mg3(PO4)2