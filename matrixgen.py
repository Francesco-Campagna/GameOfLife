from random import randint

i = int(input("Enter the dimension of i: "))
j = int(input("Enter the dimension of j: "))

with open("input.txt", "w") as file: 
    for _ in range(i):
        for _ in range(j):
            file.write(str(randint(0, 1)) + " ")
        file.write("\n")

print("Matrix generated and saved in input.txt")