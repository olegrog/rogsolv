#include <random>
#include <functional>			// for std::bind
#include <algorithm>			// for std::generate

#include "korobov.h"
#include "../containers/vel_grid.h"
#include "../workers/manager.h"

using namespace CI_griders;

const std::size_t coefficients[][1 + dimension] = {
	{    10007,  1, 1, 1,    1905,    6491,    6710,    3611,    4146,    2607,    2863 },
	{    20011,  1, 1, 1,    7191,    2057,    3758,    8928,    5960,   14809,   12988 },
	{    30011,  1, 1, 1,   13167,   26353,    2769,   26069,   14716,   14556,    8606 },
	{    50021,  1, 1, 1,   11281,    7537,   39218,   32534,   11977,    5816,   32765 }, 
	{    100003, 1, 1, 1,   20285,   68883,   49739,   25348,   68757,   93907,   46351 },
	{    200003, 1, 1, 1,   47369,  188507,   54145,  156036,  158419,   37051,   42494 },
	{    300017, 1, 1, 1,   81575,  103565,  136172,  101475,   54078,  262899,  170731 },
	{    400009, 1, 1, 1,  193141,  206577,  390670,  296791,   20804,   14959,  331221 },
	{    500009, 1, 1, 1,   42535,  193663,  307439,  182488,  487373,   37415,  418387 },
	{    750019, 1, 1, 1,   10525,  522832,  667416,  625465,  102362,  332766,  523439 },
	{   1000003, 1, 1, 1,  417564,  171019,  163483,  410620,  615303,  611111,  188073 },
	{   1500007, 1, 1, 1,  413996,  388189,  943278, 1496508,  434758,  733031,  985685 },
	{   2000003, 1, 1, 1,  832685, 1269182, 1228431,  532894,  174792,  458201,  527381 },
	{   3000017, 1, 1, 1,  368334,  166765, 2871452,  407635,  979274, 1865572, 2703215 },
	{   4000037, 1, 1, 1,   72362,  210611,   92212,  583028,  681897, 2974319, 1680656 },
	{   6000011, 1, 1, 1, 1323844, 1723313, 1392620,  251332, 5750225,  908859, 5328166 }, 
	{   8000009, 1, 1, 1,   93973, 6914802, 3957321,  907968, 4380879, 3879127, 4791477 }, 
	{  10000019, 1, 1, 1, 1833663, 3609180, 7252140, 5522715, 8914182, 4652083, 6262402 },
	{1000000000, 1, 1, 1,       1,       1,       1,       1,       1,       1,       1 }
};

inline real frac (real x)
{
	real int_part;
	return std::modf (x, &int_part);
} 

Point&& Korobov_iterator::operator* () const
{
	Point point;
	for (int i = 0; i < dimension; ++i) { 
		point[i] = frac (
			shift[i] + static_cast<real> (coefficients[line][i+1] * (iter+1)) 
				/ static_cast<real> (coefficients[line][0])
		);
	}
	return std::move (point);
}

Korobov_grid::Korobov_grid (std::size_t size)
{
	for (line = 0; coefficients[line][0] < size; ++line);
	size_ = coefficients[line][0];
	std::default_random_engine generator (1000);
	std::uniform_real_distribution<real> distrib (0, 1);
	std::generate (shift.begin (), shift.end (), std::bind (distrib, generator));
}

void Korobov::generate (real time_step)
{
	int R = mapper ().radius ();
	double vel_step = Vel_grid::cut_vel () / R;
	ci::gen (time_step, R, R, mapper (), mapper (), vel_step,
				mass, mass, particle, particle, grid);
}

void Korobov::info ()
{
	manager ().get_printer ().var ("Type", "Korobov");
	manager ().get_printer ().var ("Number of points", grid.size ());
}
