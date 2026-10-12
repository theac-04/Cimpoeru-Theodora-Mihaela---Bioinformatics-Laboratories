def count_seq(s: str) -> None:
    c_a = 0
    c_c = 0
    c_t = 0
    c_g = 0

    for char in s:
        if char == "a":
            c_a +=1
        elif char == "c":
            c_c +=1
        elif char == "t":
            c_t +=1
        elif char == "g":
            c_g +=1
        else:
            print("Wrong input")

    print(f"A= {c_a}")
    print(f"C= {c_c}")
    print(f"T= {c_t}")
    print(f"G= {c_g}")


if __name__ == "__main__":
    s = "aaaaatttttgggggccccccc"
    count_seq(s)