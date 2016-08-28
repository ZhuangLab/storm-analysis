/*
 * GPU code for homotopy calculations.
 *
 * Hazen 1/13
 */

__kernel void computeCVec(const __global float* g_matrix,
	                  const __global float* a_y_vec,
	                  const __global float* x_vec,
                          const __global int* on_indices,
	                  uint ncols,
                          uint on_indices_size,
                          __global float* c_vec)
{
    uint y = get_global_id(0);

    if (y < ncols){
        const __global float* g_row = g_matrix + y * ncols;

        int index;
        float sum = 0.0;
        for (int i = 0; i < on_indices_size; i++){
            index = on_indices[i];
            sum += g_row[index] * x_vec[index];
        }

	c_vec[y] = a_y_vec[y] - sum;
    }
}

__kernel void computeWork1(const __global float* a_matrix,
                           const __global float* d_vec,
                           const __global int* on_indices,
                           uint nrows,
			   uint ncols,
			   uint on_indices_size,
			   __global float* work1)
{
    uint y = get_global_id(0);

    if (y < nrows){
        const __global float* a_row = a_matrix + y * ncols;

        float sum = 0.0;
        for (int i = 0; i < on_indices_size; i++){
            sum += a_row[on_indices[i]]*d_vec[i];
        }

        work1[y] = sum;
    }
}

__kernel void computeG1Vec(const __global float* a_matrix_trans,
	                   const __global float* c_vec,
	                   const __global float* work1,
                           uint nrows,
			   uint ncols,
                           float lambda,
                           __global float* g1_vec)
{
    uint y = get_global_id(0);

    if (y < ncols){

        const __global float* a_row = a_matrix_trans + y * nrows;

        float sum = 0.0;
        for (int i = 0; i < nrows; i++){
            sum += a_row[i] * work1[i];
        }
	    
        g1_vec[y] = (lambda - c_vec[y])/(1.0 - sum);
    }
}

/*
 * The MIT License
 *
 * Copyright (c) 2013 Zhuang Lab, Harvard University
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
