/**
 * @file main.c
 * @author  Erez Hecht <erez.hecht@mail.huji.ac.il>
 * @version 1.0
 * @date 15 Nov 2018
 *
 * @brief Program to compare strings
 *
 * @section LICENSE
 * This program is not a free software;
 *
 * @section DESCRIPTION
 * The program analyzes proteins.
 * Input  : Proteins to analyze in pdb files.
 * Process: Read each file and produce CG, RG and Dmax values for each.
 * Output : Print to the screen for each file.
 */

//-------------------------------Includes-------------------------------
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
//-------------------------------Constants-------------------------------

/**
 * @def FILE_OPEN_ERROR "Error opening file: %s\n"
 * @brief An error message for file open errors
 */
#define FILE_OPEN_ERROR "Error opening file: %s\n"
/**
 * @def MAX_SEQUENCES 100
 * @brief The max number of sequences in a file.
 * */
#define MAX_SEQUENCES 100
/**
 * @def MAX_LINE_LEN 100
 * @brief The max number of chars in a line.
 * */
#define MAX_LINE_LEN 100

#define M_LOC 2
#define S_LOC 3
#define G_LOC 4
#define NUM_OF_ARGS 5
#define MALLOC_ERROR "Error allocating memory."
#define REALLOC_ERROR "Error re-allocating memory."
#define MAX(x, y) ((x)>(y)?(x):(y))


/**
 * TODO:
 * 1. Read the file into an array of sequences. Each sequence is preceeded by a header line
 * starting with ">". A sequence can be several lines, in which case you need to turn it into a continuous string.
 * You can assume the file has less than 100 sequences, and also no sentence is longer than 100 chars.
 * 2. Matching he sequences: 3 options for matching x1..xi, y1..yj:
 *      a) Match x1..x(i-1) to y1..y(j-1) and then:
 *          F(i,j) = F(i-1, j-1) + m : if xi==yj
 *          F(i,j) = F(i-1, j-1) + s : if xi!=yj
 *      b) Match x1..xi to y1..y(j-1) and then match between gap and yj:
 *          F(i,j) = F(i, j-1) + g
 *      c) Match x1..x(i-1) to y1..yj and then match between gap and xi:
 *          F(i,j) = F(i-1, j) + g
 * We'll chose the maximal value for F(i,j), and put it in a table.
 * Finally we'll choose F(n,m).
 * In each cell of the table well save a pointer to the previous cell that gave the best score, so that we can
 * recreate the best match.
 *
 * */
//-------------------------------Structs-------------------------------

typedef struct
{
    char *name;
    char *sequence;
} Sequence;

Sequence *initializeSequence()
{
    Sequence *seq = (Sequence *) malloc(sizeof(Sequence));
    if (seq == NULL)
    {
        fprintf(stderr, MALLOC_ERROR);
        exit(EXIT_FAILURE);
    }
    seq->name = (char *) malloc(sizeof(char));
    seq->sequence = (char *) malloc(sizeof(char));
    seq->sequence[0] = '\0';
    return seq;
}

void freeSequence(Sequence sequence)
{
    free(sequence.sequence);
    free(sequence.name);
//    free(sequence);
}

typedef struct TableCell
{
    size_t score;
    struct TableCell *previous;
} TableCell;

TableCell *initTableCell()
{
    TableCell *tc = (TableCell *) malloc(sizeof(TableCell));
    if (tc == NULL)
    {
        fprintf(stderr, "Error Mallocating");
        exit(EXIT_FAILURE);
    }
    tc->score = 0;
    tc->previous = NULL;
    return tc;
}

void freeCell(TableCell tableCell)
{
    if (tableCell.previous != NULL)
    {
        free(tableCell.previous);
    }
//    free(tableCell);
}

//-------------------------------Globals-------------------------------
Sequence sequences[MAX_SEQUENCES];
size_t numberOfSequences = 0;

/**
 * @brief Free the memory of the sequences.
 *
 * @param numberOfSequences the number of sequences to free.
 * */
void freeSequences()
{
    for (int i = 0; i < numberOfSequences; ++i)
    {
        freeSequence(sequences[i]);
    }
}
//-------------------------------Utility Functions-------------------------------

void *allocateMemory(size_t num, size_t size)
{
    void *out = malloc(num * size);
    if (out == NULL)
    {
        freeSequences();
        fprintf(stderr, MALLOC_ERROR);
        exit(EXIT_FAILURE);
    }
    return out;
}

void *reallocateMemory(void *src, size_t num, size_t size)
{
    void *out = realloc(src, num * size);
    if (out == NULL)
    {
        freeSequences();
        fprintf(stderr, REALLOC_ERROR);
        exit(EXIT_FAILURE);
    }
    return out;
}

//-------------------------------Read File-------------------------------//
/**
 * @brief Get rid of unwanted characters in the given line, and return the new size of the line.
 *
 * @param line - the line to change
 * @param lineLen - the length of the line
 * */
size_t getRidOfNewline(char *line, size_t lineLen)
{
    size_t newLen = lineLen;
    if (lineLen > 0 && line[lineLen - 1] == '\n')
    {
        line[lineLen - 1] = '\0';
        newLen--;
    }
    if (lineLen > 1 && line[newLen - 1] == '\r')
    {
        line[lineLen - 2] = '\0';
        newLen--;
    }
    return newLen;
}

/**
 * @brief sets the name parameter of the given Sequence.
 * */
void setName(Sequence *seq, char *line, size_t lineLen)
{
//    seq->name = (char *) realloc(seq->name, lineLen * sizeof(char));
    seq->name = (char *) reallocateMemory(seq->name, lineLen, sizeof(char));
    //-1 for the ">" at the start.
    strncpy(seq->name, &line[1], lineLen - 1);
}

size_t readFile(const char *address)
{
    FILE *f = fopen(address, "r");
    if (f == NULL)
    {
        fprintf(stderr, FILE_OPEN_ERROR, address);
        exit(EXIT_FAILURE);
    }
    char line[MAX_LINE_LEN];
    size_t numOfSequences = 0;
    Sequence *curr = NULL;
    while (fgets(line, MAX_LINE_LEN, f) && numOfSequences < MAX_SEQUENCES)
    {
        //Get the length of the line after we get rid of \r and \n in it:
        size_t lineLen = getRidOfNewline(line, strlen(line));
        //Case: Line is the start of a new Sequence segment:
        if (line[0] == '>')
        {
            curr = initializeSequence();
            sequences[numOfSequences] = *curr;
            setName(curr, line, lineLen);
            numOfSequences++;
        } else
        {
//            curr->sequence = (char *) realloc(curr->sequence, (strlen(curr->sequence) + lineLen) * sizeof(char) + 1);
            curr->sequence = (char *) reallocateMemory(curr->sequence, (strlen(curr->sequence) + lineLen) + 1,
                                                       sizeof(char));

            strcat(curr->sequence, line);
        }
    }
    return numOfSequences;
}

//-------------------------------Compare Sequences-------------------------------//
/**
 * TODO:
 * 1. For each two Sequences, create a 2d Table of TableCells
 * 2. Iterate over the table from top left to bottom right, and fill each spot accordingly:
 * If we are currently on table[i][j], take :
 * max{table[i-1][j], table[i][j-1], table[i-1][j-1]} and insert the value in, and furthermore save the pointer to
 * the max TableCell.
 * 3. Traverse back from Table[n][m] to get the actual best match.
 * */


TableCell *getTableCell(TableCell *table, const size_t row, const size_t col, const size_t cols)
{
    return table + (row * cols + col);
}

/**
 * @brief initializes the first row and column of a given table.
 * */
void initTable(TableCell *table, const size_t rows, const size_t cols, const int g)
{
    //Initialize the first column:
    for (size_t i = 0; i < rows; ++i)
    {
        getTableCell(table, i, 0, cols)->score = i * g;
    }
    //Initialize the first row:
    for (size_t j = 0; j < cols; ++j)
    {
        getTableCell(table, 0, j, cols)->score = j * g;
    }

}

/**
 * @brief Compares two sequences to find the longest match between them.
 *
 * @param seq1 the first sequence to compare
 * @param seq2 the second sequence to compare
 * @param m the weight of a match
 * @param s the weight of a miss-match
 * @param g the weight of a gap
 *
 * @return the length of the longest match
 * */
void compareSequences(const Sequence *seq1, const Sequence *seq2, const int m, const int s, const int g)
{
    // +1 for the extra rows and cols
    size_t rows = strlen(seq1->sequence) + 1;
    size_t cols = strlen(seq2->sequence) + 1;
//    TableCell *table = (TableCell *) malloc(rows * cols * sizeof(TableCell));
    TableCell *table = (TableCell *) allocateMemory(rows * cols, sizeof(TableCell));
    initTable(table, rows, cols, g);
    //For each cell of the table, from top left to bottom right, calculate the value.
    for (size_t i = 1; i < rows; ++i)
    {
        for (size_t j = 1; j < cols; ++j)
        {
            size_t upLeft = getTableCell(table, i - 1, j - 1, cols)->score;
            upLeft += (seq1->sequence[i - 1] == seq2->sequence[j - 1]) ? m : s;
            size_t left = getTableCell(table, i, j - 1, cols)->score + g;
            size_t up = getTableCell(table, i - 1, j, cols)->score + g;
            getTableCell(table, i, j, cols)->score = MAX(MAX(upLeft, up), left);
        }
    }
    printf("%li", getTableCell(table, (rows - 1), (cols - 1), cols)->score);
    //free the memory:
    for (size_t k = 0; k < rows; ++k)
    {
        for (size_t j = 0; j < cols; ++j)
        {
//            printf("%d ", getTableCell(table, k, j, cols)->score);
            freeCell(*getTableCell(table, k, j, cols));
        }
//        printf("\n");
    }
    free(table);
}

int *getCell(int *table, const int row, const int col, const int cols)
{
    return table + (row * cols + col);
}

void compareSeq(Sequence *seq1, Sequence *seq2, int m, int s, int g)
{
    int rows = (int) strlen(seq1->sequence) + 1;
    int cols = (int) strlen(seq2->sequence) + 1;
    int *table = (int *) malloc(rows * cols * sizeof(int));
    if (table == NULL)
    {
        printf("ERRRORROROROOROR");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < rows; ++i)
    {
        *getCell(table, i, 0, cols) = i * g;
    }
    for (int i = 0; i < cols; ++i)
    {
        *getCell(table, 0, i, cols) = i * g;
    }
    for (int i = 1; i < rows; ++i)
    {
        for (int j = 1; j < cols; ++j)
        {
            int upleft = *getCell(table, (i - 1), (j - 1), cols);
            int left = *getCell(table, i, (j - 1), cols) + g;
            int up = *getCell(table, (i - 1), j, cols) + g;
            upleft += (seq1->sequence[i - 1] == seq2->sequence[j - 1]) ? m : s;
            *getCell(table, i, j, cols) = MAX(MAX(upleft, left), up);
        }
    }
    for (int k = 0; k < rows; ++k)
    {
        for (int i = 0; i < cols; ++i)
        {
            printf("%i ", *getCell(table, (k), (i), cols));
        }
        printf("\n");
    }
    printf("Score for alignment of sequence %s to sequence %s is %i\n", seq1->name, seq2->name,
           *getCell(table, (rows - 1), (cols - 1), cols));
    free(table);
}


//V3: WORKS!!!!
int compareSeq2(Sequence *seq1, Sequence *seq2, int m, int s, int g)
{
    unsigned long int rows = strlen(seq1->sequence) + 1;
    unsigned long int cols = strlen(seq2->sequence) + 1;
    int **table = (int **) malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++)
    {
        table[i] = (int *) malloc(cols * sizeof(int));
    }
    if (table == NULL)
    {
        printf("ERRRORROROROOROR");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < rows; ++i)
    {
        table[i][0] = i * g;
    }
    for (int i = 0; i < cols; ++i)
    {
        table[0][i] = i * g;
    }
    for (int i = 1; i < rows; ++i)
    {
        for (int j = 1; j < cols; ++j)
        {
            int upleft = table[i - 1][j - 1];
            int left = table[i][j - 1] + g;
            int up = table[i - 1][j] + g;
            upleft += (seq1->sequence[i - 1] == seq2->sequence[j - 1]) ? m : s;
            table[i][j] = (int) fmax(fmax(upleft, left), up);
        }
    }
    for (int k = 0; k < rows; ++k)
    {
        for (int i = 0; i < cols; ++i)
        {
            printf("%d ", table[k][i]);
        }
        printf("\n");
    }
    return table[rows - 1][cols - 1];
}

//-------------------------------Main-------------------------------//
void validateArgs(int argc)
{
    if (argc != NUM_OF_ARGS)
    {
        printf("Usage: text.txt m s g");
        exit(EXIT_FAILURE);
    }
}

int getWeight(const char *src)
{
    char *temp = NULL;
    int weight = (int) strtol(src, &temp, 10);
    if (strlen(temp))
    {
        printf("Not a number");
        exit(EXIT_FAILURE);
    }
    return weight;

}

int main(int argc, char *argv[])
{
    validateArgs(argc);

    numberOfSequences = readFile(argv[1]);
    for (int k = 0; k < numberOfSequences; ++k)
    {
        printf("%s\n", sequences[k].sequence);
    }
    int m = getWeight(argv[M_LOC]);
    int s = getWeight(argv[S_LOC]);
    int g = getWeight(argv[G_LOC]);
    for (int i = 0; i < numberOfSequences; ++i)
    {
        for (int j = i + 1; j < numberOfSequences; ++j)
        {
//            compareSequences(&sequences[i], &sequences[j], m, s, g);
//            printf("\n-----------------------------------------------------------------\n");
            compareSeq(&sequences[i], &sequences[j], m, s, g);
            printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
//            printf("%d\n", compareSeq2(&sequences[i], &sequences[j], m, s, g));
//            printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");


        }
    }
    freeSequences();


}