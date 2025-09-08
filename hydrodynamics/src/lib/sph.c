#include <math.h>

#ifdef SPHAXES

struct particle
{
	double mass;
	double density;
	double pressure;
	double position[SPHAXES];
	double velocity[SPHAXES];
};

double sdensity(
	struct particle **particles,
	int count,
	double point[SPHAXES],
	double (*weightf)(double distance)
)
{
	double density = 0.0;
	for(int p = 0; p < count; p++)
	{
		double distance = 0.0;
		for(int a = 0; a < SPHAXES; a++)
			distance += pow(point[a] - particles[p]->position[a], 2.0);
		distance = sqrt(distance);
		double weight = weightf(distance);
		density += particles[p]->mass * weight;
	}
	return density;
}

double spressure(
	struct particle **particles,
	int count,
	double point[SPHAXES],
	double (*weightf)(double distance)
)
{
	double pressure = 0.0;
	for(int p = 0; p < count; p++)
	{
		double distance = 0.0;
		for(int a = 0; a < SPHAXES; a++)
			distance += pow(point[a] - particles[p]->position[a], 2.0);
		distance = sqrt(distance);
		double weight = weightf(distance);
		pressure += (
			particles[p]->mass *
			particles[p]->pressure /
			particles[p]->density
		) * weight;
	}
	return pressure;
}

void svelocity(
	struct particle **particles,
	int count,
	double point[SPHAXES],
	double (*weightf)(double distance),
	double result[SPHAXES]
)
{
	double velocity[SPHAXES];
	for(int a = 0; a < SPHAXES; a++) velocity[a] = 0.0;
	for(int p = 0; p < count; p++)
	{
		double distance = 0.0;
		for(int a = 0; a < SPHAXES; a++)
			distance += pow(point[a] - particles[p]->position[a], 2.0);
		distance = sqrt(distance);
		double weight = weightf(distance);
		for(int a = 0; a < SPHAXES; a++)
			velocity[a] += (
				particles[p]->mass *
				particles[p]->velocity[a] /
				particles[p]->density
			) * weight;
	}
	for(int a = 0; a < SPHAXES; a++)
		result[a] = velocity[a];
}

void spressureg(
	struct particle **particles,
	int count,
	double point[SPHAXES],
	double (*weightf)(double distance),
	void (*weightgf)(double distance, double result[SPHAXES]),
	double result[SPHAXES]
)
{
	double density = sdensity(particles, count, point, weightf);
	double pressure = spressure(particles, count, point, weightf);
	double pressureg[SPHAXES];
	for(int a = 0; a < SPHAXES; a++) pressureg[a] = 0.0;
	for(int p = 0; p < count; p++)
	{
		double distance = 0.0;
		for(int a = 0; a < SPHAXES; a++)
			distance += pow(point[a] - particles[p]->position[a], 2.0);
		distance = sqrt(distance);
		double weightg[SPHAXES];
		weightgf(distance, weightg);
		for(int a = 0; a < SPHAXES; a++)
			pressureg[a] += particles[p]->mass * (
				particles[p]->pressure /
				pow(particles[p]->density, 2.0) *
				density + pressure / density
			) * weightg[a];
	}
	for(int a = 0; a < SPHAXES; a++)
		result[a] = pressureg[a];
}

void svelocityl(
	struct particle **particles,
	int count,
	double point[SPHAXES],
	double (*weightf)(double distance),
	void (*weightlf)(double distance, double result[SPHAXES]),
	double result[SPHAXES]
)
{
	double velocity[SPHAXES];
	svelocity(particles, count, point, weightf, velocity);
	double velocityl[SPHAXES];
	for(int a = 0; a < SPHAXES; a++) velocityl[a] = 0.0;
	for(int p = 0; p < count; p++)
	{
		double distance = 0.0;
		for(int a = 0; a < SPHAXES; a++)
			distance += pow(point[a] - particles[p]->position[a], 2.0);
		distance = sqrt(distance);
		double weightl[SPHAXES];
		weightlf(distance, weightl);
		for(int a = 0; a < SPHAXES; a++)
			velocityl[a] += particles[p]->mass * (
				(particles[p]->velocity[a] - velocity[a]) /
				particles[p]->density
			) * weightl[a];
	}
	for(int a = 0; a < SPHAXES; a++)
		result[a] = velocityl[a];
}

#endif
