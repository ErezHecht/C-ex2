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
        fprintf(stderr, "Error Mallocating");
        exit(EXIT_FAILURE);
    }
    seq->name = (char *) malloc(sizeof(char));
    seq->sequence = (char *) malloc(sizeof(char));
    return seq;
}


typedef struct TableCell
{
    int score;
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
//-------------------------------Read File-------------------------------//
/**
 * @brief sets the name parameter of the given Sequence.
 * */
void setName(Sequence *seq, char *line, size_t lineLen)
{

    seq->name = (char *) realloc(seq->name, lineLen * sizeof(char));
    if (seq->name == NULL)
    {
        fprintf(stderr, "Error Mallocating");
        exit(EXIT_FAILURE);
    }

    strncpy(seq->name, &line[1], lineLen - 1);
}

unsigned int readFile(const char *address, Sequence *sequences)
{
    FILE *f = fopen(address, "r");
    if (f == NULL)
    {
        fprintf(stderr, FILE_OPEN_ERROR, address);
        exit(EXIT_FAILURE);
    }
    char line[MAX_LINE_LEN];
    unsigned int numOfSequences = 0;
    Sequence *curr = NULL;
    while (fgets(line, MAX_LINE_LEN, f) && numOfSequences < MAX_SEQUENCES)
    {
        //Next few lines are parts of sequence
        size_t lineLen = strlen(line);
        if (lineLen > 0 && line[lineLen - 1] == '\n')
        {
            line[lineLen - 1] = '\0';
            if (line[lineLen - 2] == '\r')
            {
                line[lineLen - 2] = '\0';
            }
        }
        if (line[0] == '>')
        {
            curr = initializeSequence();
            sequences[numOfSequences] = *curr;
            setName(curr, line, lineLen);
            numOfSequences++;
        } else
        {
            curr->sequence = (char *) realloc(curr->sequence, (strlen(curr->sequence) + lineLen) * sizeof(char)-1);
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

TableCell *getTableCell(TableCell *table, const int row, const int col, const int cols)
{
    return (table + (row * cols + col));
}

void initTable(TableCell *table, const int rows, const int cols, const int g)
{
    for (int i = 0; i < rows; ++i)
    {
        getTableCell(table, i, 0, cols)->score = i * g;
    }
    //Initialize the first row:
    for (int j = 0; j < cols; ++j)
    {
        getTableCell(table, 0, j, cols)->score = j * g;
    }

}

//void
//setMax(TableCell *table, const int i, const int j, const size_t cols, const int upLeft, const int up, const int left)
//{
//    TableCell *max = NULL;
//    int score = 0;
//    if (upLeft > up)
//    {
//        if (upLeft > left)
//        {
//            max = (table + (i - 1) * cols + (j - 1));
//            score = upLeft;
//        } else
//        {
//            max = (table + (i) * cols + (j - 1));
//            score = left;
//        }
//    } else
//    {
//        if (left > up)
//        {
//            max = (table + (i) * cols + (j - 1));
//            score = left;
//        } else
//        {
//            max = (table + (i - 1) * cols + (j));
//            score = up;
//        }
//    }
//    (table + i * cols + j)->previous = max;
//    (table + i * cols)->score = score;
//}
void printHelper(const Sequence *seq1, const Sequence *seq2, TableCell *firstCell, int n, int m)
{

    for (int j = -1; j <= m; ++j)
    {
        for (int i = -1; i <= n; ++i)
        {
            if (i == -1 && j == -1)
            {
                printf("   ");
            } else if (i == -1)
            {
                printf("  %c ", seq2->sequence[j - 1]);
            } else if (j == -1)
            {
                printf("  %c ", seq1->sequence[i - 1]);
            } else
            {
                printf(getTableCell(firstCell, j, i, n)->score < 0 ? "%.2d " : " %.2d ",
                       getTableCell(firstCell, j, i, n)->score);
            }

        }
        printf("\n");
    }
}

/**
 * @brief Compares two sequences to find the longest match between them.
 * */
int compareSequences(const Sequence *seq1, const Sequence *seq2, const int m, const int s, const int g)
{
    int rows = (int) strlen(seq1->sequence) + 1;
    int cols = (int) strlen(seq2->sequence) + 1;
    TableCell *table = (TableCell *) malloc(rows * cols * sizeof(TableCell));
    initTable(table, rows, cols, g);
    for (int i = 1; i < rows; ++i)
    {
        for (int j = 1; j < cols; ++j)
        {
            int upLeft = getTableCell(table, i - 1, j - 1, cols)->score;
            upLeft += (seq1->sequence[i] == seq2->sequence[j]) ? m : s;
            int left = getTableCell(table, i, j - 1, cols)->score + g;
            int up = getTableCell(table, i - 1, j, cols)->score + g;
            getTableCell(table, i, j, cols)->score = MAX(MAX(upLeft, up), left);
//            setMax(table, i, j, cols, upLeft, up, left);
        }
    }
    printHelper(seq1, seq2, table, rows, cols);
//    for (int k = 0; k < rows; ++k)
//    {
//        for (int i = 0; i < cols; ++i)
//        {
//            printf("%d ", getTableCell(table, k, i, cols)->score);
//        }
//        printf("\n");
//    }
    return getTableCell(table, (rows-1), cols, cols)->score;

}

//V2:
int compareSeq(Sequence *seq1, Sequence *seq2, int m, int s, int g)
{
    unsigned long int rows = strlen(seq1->sequence) + 1;
    unsigned long int cols = strlen(seq2->sequence) + 1;
    int *table = (int *) malloc(rows * cols * sizeof(int));
    if (table == NULL)
    {
        printf("ERRRORROROROOROR");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < rows; ++i)
    {
        *(table + (i * cols)) = i * g;
    }
    for (int i = 0; i < cols; ++i)
    {
        *(table + i) = i * g;
    }
    for (int i = 1; i < rows; ++i)
    {
        for (int j = 1; j < cols; ++j)
        {
            int upleft = *(table + (i - 1) * cols + (j - 1));
            int left = *(table + (i) * cols + (j - 1) + g);
            int up = *(table + (i - 1) * cols + j) + g;
            if (seq1->sequence[i] == seq2->sequence[j])
            {
                upleft += m;
            } else
            {
                upleft += s;
            }
            *(table + (i * cols) + j) = (int) fmax(fmax(upleft, left), up);
        }
    }
    for (int k = 0; k < rows; ++k)
    {
        for (int i = 0; i < cols; ++i)
        {
            printf("%d ", *(table + k * cols + i));
        }
        printf("\n");
    }
    return *(table + rows * cols + cols);
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
    Sequence *sequences = (Sequence *) malloc(sizeof(Sequence) * MAX_SEQUENCES);
    unsigned int numOfSequences = readFile(argv[1], sequences);
    sequences = (Sequence *) realloc(sequences, numOfSequences * sizeof(Sequence));
    for (int k = 0; k < numOfSequences; ++k)
    {
        printf("%s\n", sequences[k].sequence);
    }
    int m = getWeight(argv[M_LOC]);
    int s = getWeight(argv[S_LOC]);
    int g = getWeight(argv[G_LOC]);
    for (int i = 0; i < numOfSequences; ++i)
    {
        for (int j = i + 1; j < numOfSequences; ++j)
        {
            printf("%d\n", compareSequences(&sequences[i], &sequences[j], m, s, g));
            //printf("%d\n", compareSeq(&sequences[i], &sequences[j], m, s, g));

        }
    }


}