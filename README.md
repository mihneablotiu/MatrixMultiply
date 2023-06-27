Tema 2 ASC - Optimizare operații cu matrice - 07.05.2023
Bloțiu Mihnea-Andrei - 333CA

The general ideea of this homework was to implement the following matrix operations
C = A × B × A^t + B^t × B^t taking into consideration that A is a upper triangular matrix
and B is a random matrix.

The operation mentioned above should have been implemented in three ways:
    * One using just functions from the BLAS Atlas library
    * One using a unoptimized standard alogrithm
    * One starting from the unoptimized algorithm mentioned above and doing
    some optimizations in order to have under 12 seconds on the N = 1200 test.
and everything was implemented as follows:

1. General implementation:

    * Implementation using the BLAS Atlas library:
        - In order to implement the operations using the BLAS Atlas library I
        started with the B^T * B^T multiplication which is a operation between
        two random matrixes so I used the cblas_dgemm function for that and saved
        the result in the result_matrix;

        - After that, I made the A * B operation which is a upper triangular
        matrix multiplied by a random matrix and this can be done using the
        cblas_dtrmm function saving the result in the B matrix
        (because we no longer need it afterwards)

        - Afterwards I took the result mentioned above and multiplied it with A^T
        which is a lower triangular matrix so it can be done with the same function
        mentioned above and saved the A*B*A^T result again in the B matrix.

        - At the end, we added the result mentioned above to the result matrix
        using the cblas_daxpy function which adds two vectors in order to have
        the final result in the result_matrix and return in from the function.

        - That being said, the implementation was done using just O(N^2) additional
        space and using just functions from the BLAS Atlas libary.

    * Unoptimized Standard Alogrithm
        - In order to implement the same operation by hand, I tried to do it as
        unoptimized as possible in order to have a significant difference when doing
        the optimizations in the next algorithm;

        - That being said, each multiplication was done in a classic i-j-k way taking
        into consideration just the fact that A is a upper triangular matrix.

        - Firstly I reserved memory for each of the intermediate steps of the operation.

        - After that, to determine the A * B part, I took into consideration that it is
        a upper triangular matrix multiplied by a random matrix so the most inner loop
        should only multiply elements between the current line and the maximum number of lines;

        - Then, I had to transpose the A matrix. This was done just by changing the elements
        from i-j position to j-i position in-place for A matrix (because we do not need the
        original matrix anymore after wards) so no the A matrix becomes a lower triangular
        matrix;

        - To determine the AB * (A^T) result, I took into consideration that now (A^T) is a
        lower triangular matrix so the most inner loop should only multiply elements between
        the current column number and the maximum number of columns. Afterwards, the result
        was saved in the ABA_transpose matrix.

        - The second part of the operation was to determine the (B^T * B^T) result. Firstly,
        I transposed the B matrix (because I dont need the original one afterwards) and then
        multiply it by itself (which is a multiplication between a random matrix and itself)
        so I used the clasical i-j-k algorithm with no optimization.

        - The previous result was then added to the ABA_transpose matrix and then the final
        result was returned from the function.

        - That being said, the implementation uses 4 extra matrixes so O(4 * N^2) additional
        space with the only time optimization being the one that was necessary by the homework
        task (to take into consideration that A is a upper triangular matrix which I did by
        limiting the inner most loop to the current line/column number)

    * Optimized standard algorithm:
        - The flow of the algorithm is exactly as the unoptimized one because we were not allowed
        to change the algorithm complexity, but only to change some minor aspects to get a better
        performance out of the same flow.
        - That being said, the optimizations that I did are the following:
            - For the upper_triangular with a random matrix multiplication and for the random with
            random matrix multiplication I changed the order of the loop from i-j-k to i-k-j in
            order to get only constant and sequential access to the memory. 
                * I did not do the same thing for the random with a lower_triangular matrix multiplication
                because the most inner loop depends from the j index so j has to be before k. I did not
                find a straight-forward way to do it so unfortunatelly I could not improve this
                loop in this way.

            - For the operation mentioned above that I could not change the order of the loops in a easy
            way, I detected that result[i][j] was a constant in the most inner loop so I used a register
            variable because I wanted it to be kept in a register for a faster access to the resource

            - Also, for all the operations, where it was possible, I tried to change the direct access to
            a specific element (which implies a multiplication and a sum) to operations with only pointers
            (which implies just some incrementation) which are a lot faster than multiplications.

        - With all the optimizations mentioned above, I managed to use the same algorithm and for
        a test with 1200 * 1200 size to reduce the running time from approx. 32 seconds to 9 seconds
        which is less than 12 seconds (the value that we were asked to do in order to take the points
        for the task).

2. Cachegrind Analysis between the unoptimized and manually optimized algorithm:

    - Here I will compare the outputs of cachegrind especially for the unoptimized and optimized
    algorithm because I do not know what hides behind the BLAS functions so it would be pretty
    hard to predict why some aspects have a specific value but I will make a comparison between
    the two implementations that I know 100% what hides behind taking into considerations
    the number of Data Cache Misses/Rates, Branches Mispredicts and Instruction References:

    - First of all, taking into consideration the number of Data Cache Misses, we can see that the
    unoptimized algorithm had more than 3 times the number of data misses the optimized algorithm
    had (including read an writes). This can be explained because the optimized algorithm makes the
    multiplication of the matrixes in the i-k-j order this making only constant or sequential accesses
    to the memory so there are higher chances that the next information needed is in the cache than
    when doing a non-sequencial access. We can see that overall, the unoptimized algorithm has a 3.8%
    data miss rate while the optimized algorithm has a 3.6% miss rate.

    - Also, the information mentioned above is also valid for the Branches and Mispredicts because
    when doing only constant or sequential accesses to the memory, the spatial and temporal locality
    can take place a lot better. This would help the branch predictor to know which element is going
    to be accessed next and usually also bring the next elements to the cache memory for a faster
    future access. If the branch predictor misses because of an non-sequencial access, then we have
    a penalty because the wrong information was brought to the cache and it has to be invalidated
    and the correct information to be brought back. This aspect is also confirmed by the number of
    mispredicts that is a little bit higher at the unoptimized algorithm than at the optimized one.
    The misspredict rate is at 0.4% for the unoptimized algorithm while for the optimized algorithm
    is at 0.3%.

    - About instruction references we can see that the unoptimized algorithm has more than 3 times
    the number of instruction references the optimized algorithm has. This can be partially explained
    by the number of operations required to perform an action that is done by the two algorithms.
    
        * For example, when accessing an element in a matrix in the unoptimized version we have to do
        a multiplication (i * N) then a sum (i * N + j), then another sum (matrix + i * N + j) and then
        a dereferencing *(matrix + i * N + j)
        * When doing the same access to an element in a matrix in the optimized version we always have
        just to increment a pointer so a sum with +1 and afterwards a dereferencing.
        * This being done at each of the functions that we call in both of our algorithms, makes us
        think that it should make sense to have a lot less instruction references in the oprimized
        version than in the unoptimized version.

3. Performance comparison between the three algorithms:

    - In order to be able to compare the three algorithms I did a a few more tests than the ones
    given by the ASC team (input and input_valgrind). Those can be find in the personal_input
    file and cover a large variety of matrix sizes between 5 and 1600 (the maximum size allowed
    by the team).

    - Afterwards, I ran each of the implementations on these 21 tests and saved the results in the
    "times" directory for each of them.

    - Then, I did 2 plots in Python considering the time that took for each of the algorithms
    to compute the result for different matrix sizes. The code can be found in the plots.py
    file.

    - The plots are made by the following two aspects:
        * Firstly comparing the time necessary to compute the result by each of the algorithms
        * Secondly comparing the speedup (aka. how many times is a algorithm better than
        the others)

    * Matrix Multiplication Time Performance:
        - Here I just took the times mentioned above and plotted them;
        - We can see that all the plots have an exponential graphic but the worst being the
        unoptimized algorithm with almost 100 seconds to compute the result of a 1600x1600
        matrix multiplication;
        - After that, out manual optimizations explained above has just 23 seconds for the
        same operation which means it is probably more than 4 times better than the unoptimiez
        algorithm even if the complexity is the same.
        - The best by far is the algorithm using just BLAS functions which uses just over 2
        seconds in order to compute the same operation mentioned above.

    * Matrix Multiplication Speedup Performance:
        - Here I thought that it makes sense to look at the three implementations and make a
        plot considering how many times we managed to improve the worst algorithm. That being
        said, the only plots that make sense are:"
            * How many times is the optimized algorithm better than the unoptimized one;
            * How many times is the BLAS algorithm better than the unoptimized one;
            * How many times is the BLAS algorithm better than the optimized one;

        - To be able to plot this, in the plots.py I divided the running times of the three
        algorithms for the combinations mentioned above and we can see that:
            - The manual optimiezed algorithm is approx. 4 times better than the unoptimized
            algorithm;
            - The BLAS algorithm is approx. 10 times better than the optimiezed algorithm;
            - The BLAS algorithm is approx. more than 40 times better than the unoptimized
            algorithm;

4.  In conclusion of all the mentioned above, we did a pretty good job optimizing the standard
    algorithm but the BLAS functions are far above the optimizations I was able to think about
    in two weeks time.