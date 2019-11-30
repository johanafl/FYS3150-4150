#include "solarsystem.h"

void Solarsystem::add_planet(double mass_input, arma::vec U0)
{   
    while (array_size <= 3*num_planets + 3)
    {   // Resize all arrays if they cant fit 3 (1) more values.
        resize();
    }

    pos[3*num_planets + 0] = U0(0);
    pos[3*num_planets + 1] = U0(1);
    pos[3*num_planets + 2] = U0(2);

    vel[3*num_planets + 0] = U0(3);
    vel[3*num_planets + 1] = U0(4);
    vel[3*num_planets + 2] = U0(5);

    mass[num_planets] = mass_input;

    num_planets++;
}

arma::vec Solarsystem::get_U0()
{   
    arma::vec U0(num_planets*3*2);
    U0.zeros();

    for (int i = 0; i < num_planets*3; i++)
    {
        U0(i) = pos[i];
    }
    
    for (int i = 0; i < num_planets*3; i++)
    {
        U0(i + num_planets*3) = vel[i];
    }
    return U0;
}

void Solarsystem::resize()
{   
    array_size *= 2;
    double* tmp_pos  = new double[array_size];
    double* tmp_vel  = new double[array_size];
    double* tmp_mass = new double[array_size/3];

    for (int i = 0; i < 3*num_planets; i++)
    {
        tmp_pos[i] = pos[i];
        tmp_vel[i] = vel[i];
    }

    for (int i = 0; i < num_planets; i++)
    {
        tmp_mass[i] = mass[i];
    }

    delete[] pos;
    delete[] vel;
    delete[] mass;

    pos  = tmp_pos;
    vel  = tmp_vel;
    mass = tmp_mass;
}

arma::vec Solarsystem::acceleration(arma::vec u, double t)
{
    // double G = 1;

    arma::vec value(3);
    value.zeros();
    double r;
    double x;
    double y;
    double z;

    x = u(0);
    y = u(1);
    z = u(2);
    r = std::sqrt(x*x + y*y + z*z);

    value(0) -= 4*pi*pi*x/(r*r*r);
    value(1) -= 4*pi*pi*y/(r*r*r);
    value(2) -= 4*pi*pi*z/(r*r*r);


    return value;
}

void Solarsystem::solve_system()
{   
    int num_steps = 1e5;
    double dt = 1e-3;
    arma::vec U0 = get_U0();
    
    ForwardEuler<Solarsystem> solved(num_steps, num_planets);
    // VelocityVerlet solved(num_steps);
    solved.set_initial_conditions(U0);
    solved.solve(*this, dt);
    solved.write_to_file();   
}


// int main()
// {   
//     arma::vec U0 = {1, 0, 0, 0, 2*pi, 0};
//     double mass = 5.972e24; // Earth mass.
//     Solarsystem q;
//     q.add_planet(mass, U0);
//     q.solve_system();
//     // arma::vec initial = q.get_U0();
//     return 0;
// }



//crap:

    // void add_planet(double mass, double pos_x, double pos_y, double pos_z,
    //     double vel_x, double vel_y, double vel_z)
    // {   
    //     Planet new_planet(mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
        
    //     // resize
    //     if (num_planets == array_size)
    //     {
    //         array_size *= 2;
    //         Planet* tmp = new Planet[array_size];

    //         // double* tmp_ax = new double[array_size];
    //         // double* tmp_ay = new double[array_size];
    //         // double* tmp_az = new double[array_size];

    //         for (int i = 0; i < num_planets; i++)
    //         {
    //             tmp[i] = planets[i];
    //         }

    //         delete[] planets;
    //         // delete[] ax;
    //         // delete[] ay;
    //         // delete[] az;

    //         planets = tmp;

    //         // ax = tmp_ax;
    //         // ay = tmp_ay;
    //         // az = tmp_az;
    //     }

    //     planets[num_planets] = new_planet;
    //     num_planets += 1;
    //     acc = arma::zeros<arma::vec>(num_planets);
    // }


           // arma::vec<double> acc[];

        // for planet in planets:
        // for (int i = 0; i < num_planets; i++)
        // {
        //     for (int j = 0; j < num_planets; j++)
        //     {
        //         x = u(3*j) - u(3*i);
        //         y = u(3*j + 1) - u(3*i + 1);
        //         z = u(3*j + 2) - u(3*i + 2);
        //         r = std::sqrt(x*x + y*y + z*z);
        //         // unsure of how to get the sign correct for x, y and z.

        //         acc(3*i)     -= G*planets[j].mass*x/(r*r*r);
        //         acc(3*i + 1) -= G*planets[j].mass*y/(r*r*r);
        //         acc(3*i + 2) -= G*planets[j].mass*z/(r*r*r);

        //         // // a += G*planet.mass/|planet.r - current_planet.r|^2
        //         // // Planet MUST contain position to make this work!
        //         // r = (planet[i].r - planet[j].r);
        //         // x = (planet[i].x - planet[j].x);
        //         // y = (planet[i].y - planet[j].y);
        //         // z = (planet[i].z - planet[j].z);

        //         // ax[j] -= G*planets[i].mass*x/(r*r*r);
        //         // ay[j] -= G*planets[i].mass*y/(r*r*r);
        //         // az[j] -= G*planets[i].mass*z/(r*r*r);
        //     }
        // 