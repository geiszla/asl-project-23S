def vec_sum_flops(n):
    flops = 0

    for i in range(n - 2, -1, -1):
        flops += 6

    return flops


def vec_sum_err_flops(n):
    flops = 0

    for i in range(0, n - 1):
        flops += 6

    return flops


def vec_sum_err_branch_flops(n):
    flops = 0

    for i in range(0, n - 1):
        flops += 6

    return flops


def renormalization_flops(n, m):
    flops = vec_sum_flops(n)
    flops += vec_sum_err_branch_flops(n)

    for i in range(0, m - 1):
        flops += vec_sum_err_flops(m - i + 1)

    return flops


def multiplication_flops(k):
    flops = 3

    for n in range(1, k):
        for i in range(0, n + 1):
            flops += 3

        flops += vec_sum_flops(n**2 + n + 1)
        flops+= n**2 +2*n

    for i in range(1, k):
        flops += 3

    for i in range(0, k**2):
        flops += 2

    flops += renormalization_flops(k + 1, k)

    return flops



def multiplication(expansion_length):
    return (
        expansion_length * (expansion_length + 1) / 2 * 3
        + (expansion_length / 2) * vec_sum_flops(expansion_length)
        + (expansion_length - 1) * 2
        + expansion_length**2
        + renormalization_flops(expansion_length + 1, expansion_length)
    )


print(multiplication_flops(1547))
print(multiplication(1547))
