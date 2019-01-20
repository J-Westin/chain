#ifndef ROPE_H
#define ROPE_H




/*
namespace jw{

class Vertex {
    private:
        double x, y; // position
        double u, v; // velocity
        double a, b; // acceleration

        double rho;
        double inv_mass;

    public:
        Vertex(double x, double y, double inv_mass=1.)
        : x(x), y(x), u(0), v(0), a(0), b(0), rho(1.), inv_mass(inv_mass) {}

        double get_x()        { return x; }
        double get_y()        { return y; }
        double get_u()        { return u; }
        double get_v()        { return v; }
        double get_a()        { return a; }
        double get_b()        { return b; }
        double get_rho()      { return rho; }
        double get_inv_mass() { return inv_mass; }

        void apply_force(double f, double g){
            a += f*inv_mass;
            b += g*inv_mass;
        }

        void euler_integrate(double dt){
            x += dt*u;
            y += dt*v;
            u += dt*a;
            v += dt*b;
        }
};

class RopeJoint{
    private:
        Vertex* A, B;
        double length;

    public:
        RopeJoint(Vertex& A, Vertex& B, double length=1.)
        : A(A), B(B), length(length) {}


};

} // namespace jw
*/

#endif // ROPE_H
