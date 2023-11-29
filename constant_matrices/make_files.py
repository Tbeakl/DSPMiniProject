def output_file(i):
    with open(f"constant_{i}.txt", "w") as file:
        for _ in range(8):
            file.write(f"{i} {i} {i} {i} {i} {i} {i} {i}\n")

for i in range(1, 74):
    output_file(i)