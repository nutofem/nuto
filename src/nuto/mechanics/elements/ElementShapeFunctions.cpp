
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include <assert.h>


namespace NuTo
{
namespace ShapeFunctions1D // interval -1 to 1
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D2N(double rNaturalCoordinates, std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==2);
        rShapeFunctions[0] = 0.5*(1.-rNaturalCoordinates);
        rShapeFunctions[1] = 0.5*(1.+rNaturalCoordinates);
    }

    void DerivativeShapeFunctions1D2N(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==2);
        rDerivativeShapeFunctions[0] = -0.5;
        rDerivativeShapeFunctions[1] = 0.5;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D3N(double rNaturalCoordinates, std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==3);
        rShapeFunctions[0] = 0.5*(1.-rNaturalCoordinates)-0.5*(1.-rNaturalCoordinates*rNaturalCoordinates);
        rShapeFunctions[1] = 1.-rNaturalCoordinates*rNaturalCoordinates;
        rShapeFunctions[2] = 0.5*(1.+rNaturalCoordinates)-0.5*(1.-rNaturalCoordinates*rNaturalCoordinates);
    }

    void DerivativeShapeFunctions1D3N(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==3);
        rDerivativeShapeFunctions[0] = -0.5 + rNaturalCoordinates;
        rDerivativeShapeFunctions[1] = -2.0 * rNaturalCoordinates;
        rDerivativeShapeFunctions[2] =  0.5 + rNaturalCoordinates;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D4N(double rNaturalCoordinates, std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==4);
        double x = rNaturalCoordinates;
        rShapeFunctions[0] = -  9./16. * (x+1./3) * (x-1./3) * (x-1);
        rShapeFunctions[1] =   27./16. * (x+1.)   * (x-1./3) * (x-1);
        rShapeFunctions[2] = - 27./16  * (x+1.)   * (x+1./3) * (x-1);
        rShapeFunctions[3] =    9./16. * (x-1./3) * (x+1./3) * (x+1);

//      other option:
//        assert(rShapeFunctions.size()==4);
//        double r(rNaturalCoordinates);
//        double r2(rNaturalCoordinates*rNaturalCoordinates);
//        double r3(r2*rNaturalCoordinates);
//        rShapeFunctions[0] =  -0.0625 + 0.0625*r + 0.5625*r2 -0.5625*r3;
//        rShapeFunctions[1] =  + 0.5625 -1.6875*r -0.5625*r2 + 1.6875*r3;
//        rShapeFunctions[2] =  + 0.5625 + 1.6875*r -0.5625*r2 -1.6875*r3;
//        rShapeFunctions[3] =  -0.0625 -0.0625*r + 0.5625*r2 + 0.5625*r3;
    }

    void DerivativeShapeFunctions1D4N(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==4);
        double x = rNaturalCoordinates;
        rDerivativeShapeFunctions[0] = (-27.0*x*x + 18.0*x +  1.0) / 16.;
        rDerivativeShapeFunctions[1] = ( 81.0*x*x - 18.0*x - 27.0) / 16.;
        rDerivativeShapeFunctions[2] = (-81.0*x*x - 18.0*x + 27.0) / 16.;
        rDerivativeShapeFunctions[3] = ( 27.0*x*x + 18.0*x -  1.0) / 16.;

//      other option:
//        assert(rDerivativeShapeFunctions.size()==4);
//        double r(rNaturalCoordinates);
//        double r2(rNaturalCoordinates*rNaturalCoordinates);
//        rDerivativeShapeFunctions[0] =  + 0.0625 + 1.125*r-1.6875*r2;
//        rDerivativeShapeFunctions[1] = -1.6875-1.125*r + 5.0625*r2;
//        rDerivativeShapeFunctions[2] =  + 1.6875-1.125*r-5.0625*r2;
//        rDerivativeShapeFunctions[3] = -0.0625 + 1.125*r + 1.6875*r2;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D5N(double rNaturalCoordinates, std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==5);
        double r(rNaturalCoordinates);
        double r2(rNaturalCoordinates*rNaturalCoordinates);
        double r3(r2*rNaturalCoordinates);
        double r4(r3*rNaturalCoordinates);
        rShapeFunctions[0] =  + 0.166666666667*r -0.166666666667*r2 -0.666666666667*r3 + 0.666666666667*r4;
        rShapeFunctions[1] =  -1.33333333333*r + 2.66666666667*r2 + 1.33333333333*r3 -2.66666666667*r4;
        rShapeFunctions[2] =  + 1.0 -5.0*r2 + 4.0*r4;
        rShapeFunctions[3] =  + 1.33333333333*r + 2.66666666667*r2 -1.33333333333*r3 -2.66666666667*r4;
        rShapeFunctions[4] =  -0.166666666667*r -0.166666666667*r2 + 0.666666666667*r3 + 0.666666666667*r4;
    }

    void DerivativeShapeFunctions1D5N(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==5);
        double r(rNaturalCoordinates);
        double r2(rNaturalCoordinates*rNaturalCoordinates);
        double r3(r2*rNaturalCoordinates);
        rDerivativeShapeFunctions[0] =  + 0.166666666667-0.333333333333*r-2.0*r2 + 2.66666666667*r3;
        rDerivativeShapeFunctions[1] = -1.33333333333 + 5.33333333333*r + 4.0*r2-10.6666666667*r3;
        rDerivativeShapeFunctions[2] = -10.0*r + 16.0*r3;
        rDerivativeShapeFunctions[3] =  + 1.33333333333 + 5.33333333333*r-4.0*r2-10.6666666667*r3;
        rDerivativeShapeFunctions[4] = -0.166666666667-0.333333333333*r + 2.0*r2 + 2.66666666667*r3;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D2NSpectralOrder3(double rLocalCoordinates, std::vector<double>& rShapeFunctions)
    {
        assert(((int)rShapeFunctions.size()) == 4);
        double s2 = rLocalCoordinates*rLocalCoordinates;
        double s3 = rLocalCoordinates*s2;
        rShapeFunctions[0] = -0.125+0.125         *rLocalCoordinates      +0.625*s2-0.625*s3;
        rShapeFunctions[1] =  0.625-1.3975424859373684*rLocalCoordinates-0.625*s2+1.3975424859373684*s3;
        rShapeFunctions[2] =  0.625+1.3975424859373684*rLocalCoordinates-0.625*s2-1.3975424859373684*s3;
        rShapeFunctions[3] = -0.125-0.125         *rLocalCoordinates      +0.625*s2+0.625*s3;
    }

    void DerivativeShapeFunctions1D2NSpectralOrder3(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(((int)rDerivativeShapeFunctions.size()) == 4);
        double s2 = rLocalCoordinates*rLocalCoordinates;
        rDerivativeShapeFunctions[0] =  0.125             +1.25*rLocalCoordinates-1.875*s2;
        rDerivativeShapeFunctions[1] = -1.3975424859373684-1.25*rLocalCoordinates+4.192627457812105*s2;
        rDerivativeShapeFunctions[2] =  1.3975424859373684-1.25*rLocalCoordinates-4.192627457812105*s2;
        rDerivativeShapeFunctions[3] = -0.125             +1.25*rLocalCoordinates+1.875*s2;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D2NSpectralOrder4(double rLocalCoordinates, std::vector<double>& rShapeFunctions)
    {
        assert(((int)rShapeFunctions.size())==5);
        double s2 = rLocalCoordinates*rLocalCoordinates;
        double s3 = rLocalCoordinates*s2;
        double s4 = rLocalCoordinates*s3;
        rShapeFunctions[0] =    +0.375            *rLocalCoordinates -0.375            *s2 -0.875*s3             + 0.875*s4;
        rShapeFunctions[1] =    -1.336584577695453*rLocalCoordinates +2.041666666666666*s2 +1.336584577695453*s3 -2.041666666666666*s4;
        rShapeFunctions[2] =  1.                                      -3.333333333333333*s2                       +2.333333333333333*s4;
        rShapeFunctions[3] =    +1.336584577695453*rLocalCoordinates +2.041666666666666*s2 -1.336584577695453*s3 -2.041666666666666*s4;
        rShapeFunctions[4] =    -0.375            *rLocalCoordinates -0.375            *s2 + 0.875*s3            +0.875*s4;
    }

    void DerivativeShapeFunctions1D2NSpectralOrder4(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(((int)rDerivativeShapeFunctions.size())==5);
        double s2 = rLocalCoordinates*rLocalCoordinates;
        double s3 = rLocalCoordinates*s2;
        rDerivativeShapeFunctions[0] =  0.375             -0.75             *rLocalCoordinates-2.625               *s2+3.5*s3;
        rDerivativeShapeFunctions[1] = -1.336584577695453 +4.083333333333333*rLocalCoordinates+4.009753733086359517*s2-8.16666666666666*s3;
        rDerivativeShapeFunctions[2] =                    -6.666666666666666*rLocalCoordinates                        +9.33333333333333*s3;
        rDerivativeShapeFunctions[3] =  1.336584577695453 +4.083333333333333*rLocalCoordinates-4.009753733086359517*s2-8.16666666666666*s3;
        rDerivativeShapeFunctions[4] = -0.375             -0.75             *rLocalCoordinates+2.625               *s2+3.5*s3;
    }


}

namespace ShapeFunctions2D
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctionsPlane2D4N(const double rLocalCoordinates[2], std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==4);
        rShapeFunctions[0] = 0.25*(1.-rLocalCoordinates[0])*(1.-rLocalCoordinates[1]);
        rShapeFunctions[1] = 0.25*(1.+rLocalCoordinates[0])*(1.-rLocalCoordinates[1]);
        rShapeFunctions[2] = 0.25*(1.+rLocalCoordinates[0])*(1.+rLocalCoordinates[1]);
        rShapeFunctions[3] = 0.25*(1.-rLocalCoordinates[0])*(1.+rLocalCoordinates[1]);
    }

    void DerivativeShapeFunctionsPlane2D4N(const double rLocalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==8);
        rDerivativeShapeFunctions[0] = -0.25*(1.-rLocalCoordinates[1]);
        rDerivativeShapeFunctions[1] = -0.25*(1.-rLocalCoordinates[0]);

        rDerivativeShapeFunctions[2] = +0.25*(1.-rLocalCoordinates[1]);
        rDerivativeShapeFunctions[3] = -0.25*(1.+rLocalCoordinates[0]);

        rDerivativeShapeFunctions[4] = +0.25*(1.+rLocalCoordinates[1]);
        rDerivativeShapeFunctions[5] = +0.25*(1.+rLocalCoordinates[0]);

        rDerivativeShapeFunctions[6] = -0.25*(1.+rLocalCoordinates[1]);
        rDerivativeShapeFunctions[7] = +0.25*(1.-rLocalCoordinates[0]);
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions2D3N(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==3);
        rShapeFunctions[0] = 1.-rNaturalCoordinates[0]-rNaturalCoordinates[1];
        rShapeFunctions[1] = rNaturalCoordinates[0];
        rShapeFunctions[2] = rNaturalCoordinates[1];
    }

    void DerivativeShapeFunctions2D3N(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==6);
        rDerivativeShapeFunctions[0] = -1.0;
        rDerivativeShapeFunctions[1] = -1.0;

        rDerivativeShapeFunctions[2] = 1.0;
        rDerivativeShapeFunctions[3] = 0.0;

        rDerivativeShapeFunctions[4] = 0.0;
        rDerivativeShapeFunctions[5] = 1.0;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions2D4N(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==4);
        rShapeFunctions[0] = 0.25*(1.-rNaturalCoordinates[0])*(1.-rNaturalCoordinates[1]);
        rShapeFunctions[1] = 0.25*(1.+rNaturalCoordinates[0])*(1.-rNaturalCoordinates[1]);
        rShapeFunctions[2] = 0.25*(1.+rNaturalCoordinates[0])*(1.+rNaturalCoordinates[1]);
        rShapeFunctions[3] = 0.25*(1.-rNaturalCoordinates[0])*(1.+rNaturalCoordinates[1]);
    }

    void DerivativeShapeFunctions2D4N(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==8);
        rDerivativeShapeFunctions[0] = -0.25*(1.-rNaturalCoordinates[1]);
        rDerivativeShapeFunctions[1] = -0.25*(1.-rNaturalCoordinates[0]);

        rDerivativeShapeFunctions[2] = +0.25*(1.-rNaturalCoordinates[1]);
        rDerivativeShapeFunctions[3] = -0.25*(1.+rNaturalCoordinates[0]);

        rDerivativeShapeFunctions[4] = +0.25*(1.+rNaturalCoordinates[1]);
        rDerivativeShapeFunctions[5] = +0.25*(1.+rNaturalCoordinates[0]);

        rDerivativeShapeFunctions[6] = -0.25*(1.+rNaturalCoordinates[1]);
        rDerivativeShapeFunctions[7] = +0.25*(1.-rNaturalCoordinates[0]);
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions2D6N(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==6);
        double r(rNaturalCoordinates[0]);
        double s(rNaturalCoordinates[1]);

        rShapeFunctions[0] =  2. *(r * r + s *  s) + 4. * r * s - 3. * ( r + s ) + 1.;
        rShapeFunctions[1] =  2. * r * r - r;
        rShapeFunctions[2] =  2. * s * s -s;
        rShapeFunctions[3] = -4. * r *(r + s - 1.);
        rShapeFunctions[4] =  4. * r * s;
        rShapeFunctions[5] = -4. * s *(s + r - 1.);
    }

    void DerivativeShapeFunctions2D6N(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==12);
        double r(rNaturalCoordinates[0]);
        double s(rNaturalCoordinates[1]);

        rDerivativeShapeFunctions[0] = 4.*(r+s)-3. ;
        rDerivativeShapeFunctions[1] = 4.*(r+s)-3. ;

        rDerivativeShapeFunctions[2] = 4.*r-1.;
        rDerivativeShapeFunctions[3] = 0.;

        rDerivativeShapeFunctions[4] = 0.;
        rDerivativeShapeFunctions[5] = 4.*s-1.;

        rDerivativeShapeFunctions[6] = -8.*r-4.*s+4.;
        rDerivativeShapeFunctions[7] = -4.*r;

        rDerivativeShapeFunctions[8] = 4.*s;
        rDerivativeShapeFunctions[9] = 4.*r;

        rDerivativeShapeFunctions[10] = -4.*s;
        rDerivativeShapeFunctions[11] = -8.*s-4.*r+4.;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions2D10N(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==10);
        double r(rNaturalCoordinates[0]);
        double s(rNaturalCoordinates[1]);

        rShapeFunctions[0] =  + 1.0 -5.5*r -5.5*s + 9.0*r*r + 18.0*r*s + 9.0*s*s -4.5*r*r*r -13.5*r*r*s -13.5*r*s*s -4.5*s*s*s;
        rShapeFunctions[1] =  + 9.0*r -22.5*r*r -22.5*r*s + 13.5*r*r*r + 27.0*r*r*s + 13.5*r*s*s;
        rShapeFunctions[2] =  -4.5*r + 18.0*r*r + 4.5*r*s -13.5*r*r*r -13.5*r*r*s;
        rShapeFunctions[3] =  + 1.0*r -4.5*r*r + 4.5*r*r*r;
        rShapeFunctions[4] =  + 9.0*s -22.5*r*s -22.5*s*s + 13.5*r*r*s + 27.0*r*s*s + 13.5*s*s*s;
        rShapeFunctions[5] =  + 27.0*r*s -27.0*r*r*s -27.0*r*s*s;
        rShapeFunctions[6] =  -4.5*r*s + 13.5*r*r*s;
        rShapeFunctions[7] =  -4.5*s + 4.5*r*s + 18.0*s*s -13.5*r*s*s -13.5*s*s*s;
        rShapeFunctions[8] =  -4.5*r*s + 13.5*r*s*s;
        rShapeFunctions[9] =  + 1.0*s -4.5*s*s + 4.5*s*s*s;
    }

    void DerivativeShapeFunctions2D10N(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==20);
        double r(rNaturalCoordinates[0]);
        double s(rNaturalCoordinates[1]);

        rDerivativeShapeFunctions[0] = -5.5 + 18.0*r + 18.0*s-13.5*r*r-27.0*r*s-13.5*s*s;
        rDerivativeShapeFunctions[1] = -5.5 + 18.0*r + 18.0*s-13.5*r*r-27.0*r*s-13.5*s*s;
        rDerivativeShapeFunctions[2] =  + 9.0-45.0*r-22.5*s + 40.5*r*r + 54.0*r*s + 13.5*s*s;
        rDerivativeShapeFunctions[3] = -22.5*r + 27.0*r*r + 27.0*r*s;
        rDerivativeShapeFunctions[4] = -4.5 + 36.0*r + 4.5*s-40.5*r*r-27.0*r*s;
        rDerivativeShapeFunctions[5] =  + 4.5*r-13.5*r*r;
        rDerivativeShapeFunctions[6] =  + 1.0-9.0*r + 13.5*r*r;
        rDerivativeShapeFunctions[7] = 0.;
        rDerivativeShapeFunctions[8] = -22.5*s + 27.0*r*s + 27.0*s*s;
        rDerivativeShapeFunctions[9] =  + 9.0-22.5*r-45.0*s + 13.5*r*r + 54.0*r*s + 40.5*s*s;
        rDerivativeShapeFunctions[10] =  + 27.0*s-54.0*r*s-27.0*s*s;
        rDerivativeShapeFunctions[11] =  + 27.0*r-27.0*r*r-54.0*r*s;
        rDerivativeShapeFunctions[12] = -4.5*s + 27.0*r*s;
        rDerivativeShapeFunctions[13] = -4.5*r + 13.5*r*r;
        rDerivativeShapeFunctions[14] =  + 4.5*s-13.5*s*s;
        rDerivativeShapeFunctions[15] = -4.5 + 4.5*r + 36.0*s-27.0*r*s-40.5*s*s;
        rDerivativeShapeFunctions[16] = -4.5*s + 13.5*s*s;
        rDerivativeShapeFunctions[17] = -4.5*r + 27.0*r*s;
        rDerivativeShapeFunctions[18] = 0.;
        rDerivativeShapeFunctions[19] =  + 1.0-9.0*s + 13.5*s*s;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions2D15N(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==15);
        double r(rNaturalCoordinates[0]);
        double s(rNaturalCoordinates[1]);

        rShapeFunctions[0] =  + 1.0 -8.33333333333*r -8.33333333333*s + 23.3333333333*r*r + 46.6666666667*r*s + 23.3333333333*s*s -26.6666666667*r*r*r -80.0*r*r*s -80.0*r*s*s -26.6666666667*s*s*s + 10.6666666667*s*s*s*s + 42.6666666667*r*s*s*s + 64.0*r*r*s*s + 42.6666666667*r*r*r*s + 10.6666666667*r*r*r*r;
        rShapeFunctions[1] =  + 16.0*r -69.3333333333*r*r -69.3333333333*r*s + 96.0*r*r*r + 192.0*r*r*s + 96.0*r*s*s -42.6666666667*r*s*s*s -128.0*r*r*s*s -128.0*r*r*r*s -42.6666666667*r*r*r*r;
        rShapeFunctions[2] =  -12.0*r + 76.0*r*r + 28.0*r*s -128.0*r*r*r -144.0*r*r*s -16.0*r*s*s + 64.0*r*r*s*s + 128.0*r*r*r*s + 64.0*r*r*r*r;
        rShapeFunctions[3] =  + 5.33333333333*r -37.3333333333*r*r -5.33333333333*r*s + 74.6666666667*r*r*r + 32.0*r*r*s -42.6666666667*r*r*r*s -42.6666666667*r*r*r*r;
        rShapeFunctions[4] =  -1.0*r + 7.33333333333*r*r -16.0*r*r*r + 10.6666666667*r*r*r*r;
        rShapeFunctions[5] =  + 16.0*s -69.3333333333*r*s -69.3333333333*s*s + 96.0*r*r*s + 192.0*r*s*s + 96.0*s*s*s -42.6666666667*s*s*s*s -128.0*r*s*s*s -128.0*r*r*s*s -42.6666666667*r*r*r*s;
        rShapeFunctions[6] =  + 96.0*r*s -224.0*r*r*s -224.0*r*s*s + 128.0*r*s*s*s + 256.0*r*r*s*s + 128.0*r*r*r*s;
        rShapeFunctions[7] =  -32.0*r*s + 160.0*r*r*s + 32.0*r*s*s -128.0*r*r*s*s -128.0*r*r*r*s;
        rShapeFunctions[8] =  + 5.33333333333*r*s -32.0*r*r*s + 42.6666666667*r*r*r*s;
        rShapeFunctions[9] =  -12.0*s + 28.0*r*s + 76.0*s*s -16.0*r*r*s -144.0*r*s*s -128.0*s*s*s + 64.0*s*s*s*s + 128.0*r*s*s*s + 64.0*r*r*s*s;
        rShapeFunctions[10] =  -32.0*r*s + 32.0*r*r*s + 160.0*r*s*s -128.0*r*s*s*s -128.0*r*r*s*s;
        rShapeFunctions[11] =  + 4.0*r*s -16.0*r*r*s -16.0*r*s*s + 64.0*r*r*s*s;
        rShapeFunctions[12] =  + 5.33333333333*s -5.33333333333*r*s -37.3333333333*s*s + 32.0*r*s*s + 74.6666666667*s*s*s -42.6666666667*s*s*s*s -42.6666666667*r*s*s*s;
        rShapeFunctions[13] =  + 5.33333333333*r*s -32.0*r*s*s + 42.6666666667*r*s*s*s;
        rShapeFunctions[14] =  -1.0*s + 7.33333333333*s*s -16.0*s*s*s + 10.6666666667*s*s*s*s;
    }

    void DerivativeShapeFunctions2D15N(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==30);
        double r(rNaturalCoordinates[0]);
        double s(rNaturalCoordinates[1]);

        rDerivativeShapeFunctions[0] = -8.33333333333 + 46.6666666667*r + 46.6666666667*s-80.0*r*r-160.0*r*s-80.0*s*s + 42.6666666667*s*s*s + 128.0*r*s*s + 128.0*r*r*s + 42.6666666667*r*r*r;
        rDerivativeShapeFunctions[1] = -8.33333333333 + 46.6666666667*r + 46.6666666667*s-80.0*r*r-160.0*r*s-80.0*s*s + 42.6666666667*s*s*s + 128.0*r*s*s + 128.0*r*r*s + 42.6666666667*r*r*r;
        rDerivativeShapeFunctions[2] =  + 16.0-138.666666667*r-69.3333333333*s + 288.0*r*r + 384.0*r*s + 96.0*s*s-42.6666666667*s*s*s-256.0*r*s*s-384.0*r*r*s-170.666666667*r*r*r;
        rDerivativeShapeFunctions[3] = -69.3333333333*r + 192.0*r*r + 192.0*r*s-128.0*r*s*s-256.0*r*r*s-128.0*r*r*r;
        rDerivativeShapeFunctions[4] = -12.0 + 152.0*r + 28.0*s-384.0*r*r-288.0*r*s-16.0*s*s + 128.0*r*s*s + 384.0*r*r*s + 256.0*r*r*r;
        rDerivativeShapeFunctions[5] =  + 28.0*r-144.0*r*r-32.0*r*s + 128.0*r*r*s + 128.0*r*r*r;
        rDerivativeShapeFunctions[6] =  + 5.33333333333-74.6666666667*r-5.33333333333*s + 224.0*r*r + 64.0*r*s-128.0*r*r*s-170.666666667*r*r*r;
        rDerivativeShapeFunctions[7] = -5.33333333333*r + 32.0*r*r-42.6666666667*r*r*r;
        rDerivativeShapeFunctions[8] = -1.0 + 14.6666666667*r-48.0*r*r + 42.6666666667*r*r*r;
        rDerivativeShapeFunctions[9] = 0.;
        rDerivativeShapeFunctions[10] = -69.3333333333*s + 192.0*r*s + 192.0*s*s-128.0*s*s*s-256.0*r*s*s-128.0*r*r*s;
        rDerivativeShapeFunctions[11] =  + 16.0-69.3333333333*r-138.666666667*s + 96.0*r*r + 384.0*r*s + 288.0*s*s-170.666666667*s*s*s-384.0*r*s*s-256.0*r*r*s-42.6666666667*r*r*r;
        rDerivativeShapeFunctions[12] =  + 96.0*s-448.0*r*s-224.0*s*s + 128.0*s*s*s + 512.0*r*s*s + 384.0*r*r*s;
        rDerivativeShapeFunctions[13] =  + 96.0*r-224.0*r*r-448.0*r*s + 384.0*r*s*s + 512.0*r*r*s + 128.0*r*r*r;
        rDerivativeShapeFunctions[14] = -32.0*s + 320.0*r*s + 32.0*s*s-256.0*r*s*s-384.0*r*r*s;
        rDerivativeShapeFunctions[15] = -32.0*r + 160.0*r*r + 64.0*r*s-256.0*r*r*s-128.0*r*r*r;
        rDerivativeShapeFunctions[16] =  + 5.33333333333*s-64.0*r*s + 128.0*r*r*s;
        rDerivativeShapeFunctions[17] =  + 5.33333333333*r-32.0*r*r + 42.6666666667*r*r*r;
        rDerivativeShapeFunctions[18] =  + 28.0*s-32.0*r*s-144.0*s*s + 128.0*s*s*s + 128.0*r*s*s;
        rDerivativeShapeFunctions[19] = -12.0 + 28.0*r + 152.0*s-16.0*r*r-288.0*r*s-384.0*s*s + 256.0*s*s*s + 384.0*r*s*s + 128.0*r*r*s;
        rDerivativeShapeFunctions[20] = -32.0*s + 64.0*r*s + 160.0*s*s-128.0*s*s*s-256.0*r*s*s;
        rDerivativeShapeFunctions[21] = -32.0*r + 32.0*r*r + 320.0*r*s-384.0*r*s*s-256.0*r*r*s;
        rDerivativeShapeFunctions[22] =  + 4.0*s-32.0*r*s-16.0*s*s + 128.0*r*s*s;
        rDerivativeShapeFunctions[23] =  + 4.0*r-16.0*r*r-32.0*r*s + 128.0*r*r*s;
        rDerivativeShapeFunctions[24] = -5.33333333333*s + 32.0*s*s-42.6666666667*s*s*s;
        rDerivativeShapeFunctions[25] =  + 5.33333333333-5.33333333333*r-74.6666666667*s + 64.0*r*s + 224.0*s*s-170.666666667*s*s*s-128.0*r*s*s;
        rDerivativeShapeFunctions[26] =  + 5.33333333333*s-32.0*s*s + 42.6666666667*s*s*s;
        rDerivativeShapeFunctions[27] =  + 5.33333333333*r-64.0*r*s + 128.0*r*s*s;
        rDerivativeShapeFunctions[28] = 0.;
        rDerivativeShapeFunctions[29] = -1.0 + 14.6666666667*s-48.0*s*s + 42.6666666667*s*s*s;
    }
}

namespace ShapeFunctions3D
{

}

}
