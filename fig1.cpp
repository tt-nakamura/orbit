#include "Planet.h"
#include "const.h"
#include<cmath>
#include<fstream>

main() {
	int	i,j,y,m;
	double d;
    double v_accel(pow(GM*5e-5,2));
    double v_angle(pow(DEGREE,2));
    double x[6];
    std::vector<Orbit> orb;
    std::vector<ObsData> obs;
    std::ofstream f("fig1.txt");

	i = read(obs, "Ceres.txt");
    std::cout << i << std::endl;

    ObsData o[3];
    o[0] = obs[0];
    o[1] = obs[10];
    o[2] = obs[20];
	determine_by_3obs(orb, o);

    Planet p(orb[0], obs[0].t);
    for(i=1; i<obs.size(); i++) {
        p.update(obs[i], v_accel, v_angle);
        p.get_param(x);
        obs[i].get_date(y,m,d);
        f << y << '/' << m << '/' << d << '\t';
        for(j=1; j<3; j++) f << x[j] << '\t';
        f << x[3]/DEGREE << '\t';
        f << x[4]/DEGREE << '\t';
        f << x[5]/DEGREE << '\n';
    }
}
