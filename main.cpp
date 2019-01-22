#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <math.h>
using namespace std;

# define PI 3.14159265358979323846 

static const double w_value = 0.9; // inertia
static const double c1 = 0.3; // pbest effect
static const double c2 = 0.5; // global effect

static const int max_particle_count = 2000; // count of particle
static const int max_dimension = 20; // dimension
static const int max_iteration = 200; // iteration_time
static const int velocity_max = 800; // max and min
static const int velocity_min = -800; // max and min
static const int pos_max = 500;
static const int pos_min = -500;

double gbest_pos[max_dimension];
double gbest_cost;

typedef struct particle
{
	double pos[max_dimension];
	double velocity[max_dimension];
	double cost;

	double pbest_pos[max_dimension];
	double pbest_cost;
}Particle;

double sphere_cost(double* pos)
{
	double sum = 0.0;

	for (int i = 0; i < max_dimension; i++)
	{
		sum += pos[i] * pos[i];
	}

	return -sum;
}

double rastrigin_cost(double *pos)
{
	double sum = 0.0;

	for (int i = 0; i < max_dimension; i++)
	{
		sum += pos[i] * pos[i] - (10 * cos(2 * PI * pos[i]));
	}

	return -sum;
}

double schwefel_cost(double *pos)
{
	double sum = 0.0;

	for (int i = 0; i < max_dimension; i++)
	{
		sum -= pos[i] * sin(sqrt(abs(pos[i])));
	}

	return -sum;
}

int particle_swarm_optimization(void)
{
	string filePath = "result.txt";
	ofstream writeFile(filePath.data());

	srand((unsigned int)time(NULL));
	Particle particle_set[max_particle_count];

	//init particle
	for (int particle_count = 0; particle_count < max_particle_count; particle_count++)
	{
		for (int i = 0; i < max_dimension; i++)
		{
			particle_set[particle_count].pos[i] = ((double)rand() / RAND_MAX) * (pos_max - pos_min) + pos_min;
			particle_set[particle_count].pbest_pos[i] = particle_set[particle_count].pos[i];
			particle_set[particle_count].velocity[i] = ((double)rand() / RAND_MAX) * (velocity_max - velocity_min) + velocity_min;
		}
		particle_set[particle_count].cost = sphere_cost(particle_set[particle_count].pos);
		particle_set[particle_count].pbest_cost = particle_set[particle_count].cost;

		if (particle_count == 0)
		{
			gbest_cost = particle_set[particle_count].pbest_cost;
			for (int i = 0; i < max_dimension; i++)
				gbest_pos[i] = particle_set[particle_count].pos[i];
		}
		else
		{
			if (particle_set[particle_count].pbest_cost > gbest_cost)
			{
				gbest_cost = particle_set[particle_count].pbest_cost;
				for (int i = 0; i < max_dimension; i++)
					gbest_pos[i] = particle_set[particle_count].pos[i];
			}
		}
	}

	//best write
	double best_value = 0.0;
	// if schwefel_cost
	//double best_value = 420.9687;
	for (int i = 0; i < max_dimension - 1; i++)
		writeFile << best_value << ",";
	writeFile << best_value << "\n";

	//iter start
	for (int iter = 0; iter < max_iteration; iter++)
	{
		printf("iter : %d ", iter);
		for (int i = 0; i < max_dimension - 1; i++)
		{
			printf("%.2f, ", gbest_pos[i]);
		}
		printf("%.2f - ", gbest_pos[max_dimension - 1]);

		//this section caculate each particle cost and updata gbest & pbest
		for (int particle_count = 0; particle_count < max_particle_count; particle_count++)
		{
			particle_set[particle_count].cost = sphere_cost(particle_set[particle_count].pos);

			//pbest update
			if (particle_set[particle_count].cost > particle_set[particle_count].pbest_cost)
			{
				particle_set[particle_count].pbest_cost = particle_set[particle_count].cost;
				for (int i = 0; i < max_dimension; i++)
					particle_set[particle_count].pbest_pos[i] = particle_set[particle_count].pos[i];
			}

			//gbest update
			if (particle_set[particle_count].pbest_cost > gbest_cost)
			{
				gbest_cost = particle_set[particle_count].pbest_cost;
				for (int i = 0; i < max_dimension; i++)
					gbest_pos[i] = particle_set[particle_count].pos[i];
			}
			for (int i = 0; i < max_dimension - 1; i++)
				writeFile << particle_set[particle_count].pos[i] << ",";
			writeFile << particle_set[particle_count].pos[max_dimension - 1] << "\n";
		}

		for (int i = 0; i < max_dimension - 1; i++)
			writeFile << gbest_pos[i] << ",";
		writeFile << gbest_pos[max_dimension - 1] << "\n";
		
		////나중에 지워
		//for (int particle_count = 0; particle_count < max_particle_count; particle_count++)
		//{
		//	for (int i = 0; i < max_dimension - 1; i++)
		//		writeFile << particle_set[particle_count].pos[i] << ",";
		//	writeFile << particle_set[particle_count].pos[max_dimension - 1] << "\n";
		//}
		//for (int i = 0; i < max_dimension - 1; i++)
		//	writeFile << gbest_pos[i] << ",";
		//writeFile << gbest_pos[max_dimension - 1] << "\n";

		//this section caculate velocity and update it
		int speed_count = 0;
		for (int particle_count = 0; particle_count < max_particle_count; particle_count++)
		{
			double velocity_size = 0;

			for (int i = 0; i < max_dimension; i++)
			{
				particle_set[particle_count].velocity[i] = w_value * particle_set[particle_count].velocity[i] +
					c1 * ((double)rand() / RAND_MAX)*(particle_set[particle_count].pbest_pos[i] - particle_set[particle_count].pos[i]) +
					c2 * ((double)rand() / RAND_MAX)*(gbest_pos[i] - particle_set[particle_count].pos[i]);
				velocity_size += particle_set[particle_count].velocity[i] * particle_set[particle_count].velocity[i];
			}

			velocity_size = sqrt(velocity_size);
			if (velocity_size > velocity_max)
			{
				speed_count++;
				double temp_value = velocity_max / velocity_size;
				for (int i = 0; i < max_dimension; i++)
					particle_set[particle_count].velocity[i] *= temp_value;
			}
			

			//position update
			for (int i = 0; i < max_dimension; i++)
			{
				particle_set[particle_count].pos[i] += particle_set[particle_count].velocity[i];
				if (particle_set[particle_count].pos[i] > pos_max)
					particle_set[particle_count].pos[i] = pos_max;
				else if (particle_set[particle_count].pos[i] < pos_min)
					particle_set[particle_count].pos[i] = pos_min;
			}
		}
		double percent = (double)speed_count / (double)max_particle_count;
		printf("max count : %.2f %% \n", percent*100);
	}

	writeFile.close();
	return 0;
}

int main(void)
{
	particle_swarm_optimization();
	system("pause");

	return 0;
}