num_files = 9 # Number of output files to read
first_line = 18 # First line of the output file to read

for i in range(num_files):
    # Empty the output file
    op_f = open(f"outputCleaned{i+1}.txt", "w")
    op_f.close()

for i in range(num_files):
    f = open(f"output{i+1}.txt", "r")
    lines = f.readlines()

    count = 1
    for line in lines:
        if count < first_line:
            count += 1
            continue
        
        op_f = open(f"outputCleaned{i+1}.txt", "a")

        # Write the first, second, last and second last columns
        op_f.write(f"{line.split()[0]} {line.split()[1]} {line.split()[-1]} {line.split()[-2]}\n")
        op_f.close()
        count += 1
        