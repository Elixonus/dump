#pragma once

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
);

double spressure(
	struct particle **particles,
	int count,
	double point[SPHAXES],
	double (*weightf)(double distance)
);

void svelocity(
	struct particle **particles,
	int count,
	double point[SPHAXES],
	double (*weightf)(double distance),
	double result[SPHAXES]
);

void spressureg(
	struct particle **particles,
	int count,
	double point[SPHAXES],
	double (*weightf)(double distance),
	void (*weightgf)(double distance, double result[SPHAXES]),
	double result[SPHAXES]
);

void svelocityl(
	struct particle **particles,
	int count,
	double point[SPHAXES],
	double (*weightf)(double distance),
	void (*weightlf)(double distance, double result[SPHAXES]),
	double result[SPHAXES]
);

#endif
