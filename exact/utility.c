/* Utility functions for modeling */
#ifndef __UTILITY_FUNCTIONS__
#define __UTILITY_FUNCTIONS__

#include "mex.h"
#include "matrix.h"
#include <math.h>

/* Very simple log1p function since Windows doesn't have one */
/* Taken and modified from http://www.johndcook.com/cpp_log_one_plus_x.html */
#ifndef __GNUC__
inline double log1p(double x)
{
	assert(x > -1.0);

    if (x > 1e-4 || x < -1e-4)
    {
        /* x is large enough that the obvious evaluation is OK */
        return log(1.0 + x);
    }

    /* Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
       Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8 */

    return (-0.5 * x + 1.0) * x;
}
#endif

/* Generic comparator for double */
int doubleCompare(const void *arg1, const void *arg2)
{
	double* num1 = (double*)arg1;
	double* num2 = (double*)arg2;
	double difference = *num1 - *num2;
	if (difference > 0)
	{
		return 1;
	}
	else if (difference < 0)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

/* Returns index into 2D array; matches MATLAB convention */
inline size_t index2d(size_t d1, size_t d1Size, size_t d2)
{
	assert(d1 < d1Size);
	return d1 + d2 * d1Size;
}

/* Returns index into 3D array; matches MATLAB convention */
inline size_t index3d(size_t d1, size_t d1Size, size_t d2, size_t d2Size, size_t d3)
{
	assert(d1 < d1Size && d2 < d2Size);
	return d1 + d2 * d1Size + d3 * d1Size * d2Size;
}

/* Returns the maximum int in the given data array */
inline unsigned int findmax(unsigned int* data, size_t size)
{
	unsigned int i, value = 0;
	for (i = 0; i < size; ++i)
	{
		/* Find maxNumSents */
		if (data[i] > value)
		{
			value = data[i];
		}
	}
	return value;
}

/* Sums the given data array */
inline unsigned int sum(unsigned int* data, size_t size)
{
	unsigned int i, value = 0;
	for (i = 0; i < size; ++i)
	{
		value += data[i];
	}
	return value;
}

/* Samples multiple times from normalized multinomial using vector of log probs */
/* Note: Destroys rand; buffer should be of same size as logProb, and is used as scratch space; if buffer is NULL, logProb is destroyed */
/* TODO - temperature implementation seems glitchy because it is not invariant wrt multiplicative scaling of log probs; for now do not use */
void sample_multiple_normalized_log(unsigned int* samples, unsigned int numSamples, double* logProb, double* buffer, size_t count, double inverseTemperature, double* rand)
{
	double sum = 0.0;
	unsigned int i;
	unsigned int sampleIdx;
	double maxLogProb = -DBL_MAX;
	unsigned int maxLogProbPosition = 0;

	/* HACK HACK - temporary measure since temperature scaling is wrong */
	assert(inverseTemperature == 1.0);

	if (buffer == NULL)
	{
		buffer = logProb;
	}
	for (i = 0; i < count; ++i)
	{
		buffer[i] = logProb[i] * (inverseTemperature != DBL_MAX ? inverseTemperature : 1.0);
		if (buffer[i] > maxLogProb)
		{
			maxLogProbPosition = i;
			maxLogProb = buffer[i];
		}
	}
	if (inverseTemperature == DBL_MAX)
	{
		for (i = 0; i < numSamples; ++i)
		{
			samples[i] = maxLogProbPosition;
		}
		return;
	}
	for (i = 0; i < count; ++i)
	{
		buffer[i] = exp(buffer[i] - maxLogProb);
		sum += buffer[i];
	}

	/* Rescale random numbers to fit, and sort them */
	for (i = 0; i < numSamples; ++i)
	{
		assert(rand[i] >= 0 && rand[i] <= 1);
		rand[i] *= sum;
	}
	qsort(rand, numSamples, sizeof(double), doubleCompare);

	sum = 0.0;
	sampleIdx = 0;
	for (i = 0; i < count; ++i)
	{
		sum += buffer[i];
		for (; sampleIdx < numSamples && rand[sampleIdx] < sum; ++sampleIdx)
		{
			samples[sampleIdx] = i;
		}
		if (sampleIdx == numSamples)
		{
			break;
		}
	}

	assert(sampleIdx == numSamples);
}

/* Samples from normalized multinomial using vector of log probs */
/* Note: destroys logProb if buffer is NULL */
unsigned int sample_normalized_log(double* logProb, double* buffer, size_t count, double inverseTemperature, double rand)
{
	unsigned int sample;
	sample_multiple_normalized_log(&sample, 1, logProb, buffer, count, inverseTemperature, &rand);
	return sample;
}

/* Samples multiple times from normalized multinomial using vector of regular probs */
/* Destroys rand */
void sample_multiple_normalized(unsigned int* samples, unsigned int numSamples, double* prob, size_t count, double inverseTemperature, double* rand)
{
	double sum = 0.0;
	unsigned int i;
	unsigned int sampleIdx;

	/* HACK HACK - temporary measure since temperature scaling is wrong */
	assert(inverseTemperature == 1.0);

	for (i = 0; i < count; ++i)
	{
		sum += prob[i];
	}

	/* Rescale random numbers to fit, and sort them */
	for (i = 0; i < numSamples; ++i)
	{
		assert(rand[i] >= 0 && rand[i] <= 1);
		rand[i] *= sum;
	}
	qsort(rand, numSamples, sizeof(double), doubleCompare);

	sum = 0.0;
	sampleIdx = 0;
	for (i = 0; i < count; ++i)
	{
		sum += prob[i];
		for (; sampleIdx < numSamples && rand[sampleIdx] < sum; ++sampleIdx)
		{
			samples[sampleIdx] = i;
		}
		if (sampleIdx == numSamples)
		{
			break;
		}
	}

	assert(sampleIdx == numSamples);
}

/* Samples from normalized multinomial using vector of probs */
unsigned int sample_normalized(double* prob, size_t count, double inverseTemperature, double rand)
{
	unsigned int sample;
	sample_multiple_normalized(&sample, 1, prob, count, inverseTemperature, &rand);
	return sample;
}

/* Computes \log(\sum_i \exp(x_i)) in a way that reduces underflow */
inline double sum_log_probs(double* x, size_t count)
{
	/* Find max log prob */
	double max = x[0];
	unsigned int maxLocation = 0;
	unsigned int i;
	double sum = 0.0;
	for (i = 1; i < count; ++i)
	{
		if (max > x[i])
		{
			max = x[i];
			maxLocation = i;
		}
	}

	/* Perform sum of all the terms that are not the max, then log1p to add the max */
	for (i = 0; i < count; ++i)
	{
		if (i != maxLocation)
		{
			sum += exp(x[i] - max);
		}
	}

	return log1p(sum) + max;
}

/* Computes \log(\exp(x) + \exp(y)) in a way that reduces underflow */
inline double sum_two_log_probs(double x, double y)
{
	if (x > y)
	{
		return log1p(exp(y - x)) + x;
	}
	else
	{
		return log1p(exp(x - y)) + y;
	}
}

inline void drawnow(void)
{
	mexEvalString("drawnow;");
}

inline void tic(void)
{
	mexCallMATLAB(0, NULL, 0, NULL, "tic");
}

inline void toc(void)
{
	mexCallMATLAB(0, NULL, 0, NULL, "toc");
	drawnow();
}

inline void fast_print(char* str)
{
	mexPrintf(str);
	drawnow();
}

inline void log_start(char* name)
{
	tic();
	mexPrintf("Entering %s\n", name);
	drawnow();
}

inline void log_end(char* name)
{
	toc();
	mexPrintf("Leaving %s\n", name);
	drawnow();
}

#endif
