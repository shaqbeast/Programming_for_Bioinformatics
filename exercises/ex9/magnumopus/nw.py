def create_matrix(seq_a: str, seq_b: str):
    matrix_score = []
    for i in range(len(seq_a) + 1):
        matrix_score.append(list())
        for _ in range(len(seq_b) + 1):
            matrix_score[i].append(0)
    
    return matrix_score

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    matrix_score = create_matrix(seq_b, seq_a) # score matrix 

    # directions matrix 
    # T = top
    # L = left
    # D = diagonal
    matrix_pointers = create_matrix(seq_b, seq_a)
    matrix_pointers[0][0] = None

    # populating row 1 for both matrices
    initial_score = 0
    p_index = 1 # pointers index
    for column in range(len(matrix_score[0])):
        matrix_score[0][column] = initial_score
        if column + 1 != len(matrix_score[0]):
            matrix_pointers[0][column + 1] = "L"
        initial_score += -1

    # populating column 1 for both matrices 
    initial_score = 0
    p_index = 1
    for row in range(len(matrix_score)):
        matrix_score[row][0] = initial_score
        if row + 1 != len(matrix_score):
            matrix_pointers[row + 1][0] = "T"
        initial_score += -1 

    # calculating scores for entire matrix
    # i = row
    # j = column
    for i in range(1, len(matrix_score)):
        for j in range(1, len(matrix_score[i])):
            if seq_b[i - 1] == seq_a[j - 1]: # if both sequences at that position have the same nucleotide
                diagonal = matrix_score[i - 1][j - 1] + match 
            else:
                diagonal = matrix_score[i - 1][j - 1] + mismatch 
            left = matrix_score[i][j - 1] + gap 
            top = matrix_score[i - 1][j] + gap 
            matrix_score[i][j] = max(diagonal, left, top) # calculates max value

            # calculating pointers matrix
            matrix_pointers[i][j] = []
            if matrix_score[i][j] == diagonal:
                matrix_pointers[i][j].append("D")
            if matrix_score[i][j] == left:
                matrix_pointers[i][j].append("L")
            if matrix_score[i][j] == top:
                matrix_pointers[i][j].append("T")
            
    # create sequence with values from pointers matrix
    # backtracing portion 
    n = len(matrix_score) - 1
    m = len(matrix_score[0]) - 1
    current_pointer = matrix_pointers[n][m] # bottomright corner cell 
    score = matrix_score[n][m]
    horizontal = ""
    vertical = ""
    index_h = len(seq_a) - 1
    index_v = len(seq_b) - 1
    while current_pointer != None:
        if current_pointer[0] == 'D':
            n -= 1
            m -= 1
            current_pointer = matrix_pointers[n][m] # move
            horizontal += seq_a[index_h]
            vertical += seq_b[index_v]
            index_h -= 1
            index_v -= 1
        elif current_pointer[0] == 'L':
            m -= 1
            current_pointer = matrix_pointers[n][m] # move
            horizontal += seq_a[index_h]
            vertical += '-'
            index_h -= 1
        else:
            n -= 1
            current_pointer = matrix_pointers[n][m] # move
            horizontal += '-'
            vertical += seq_b[index_v]
            index_v -= 1

    # reverse the sequences 
    aligned_1 = horizontal[::-1]
    aligned_2 = vertical[::-1]

    
    
    return (aligned_1, aligned_2), score
