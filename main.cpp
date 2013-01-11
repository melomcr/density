/*
    This program is a basis for comparison of different functions that
    generate random numbers according to Tsallis distribution.
    Copyright (C) 2013 Marcelo C. R. Melo

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <fstream>
#include <list>
#include <exception>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/random.hpp>

////// Structures

struct SimData {
    
    int NumLoops;
    float qa, qt, qv;
    
    int myseed ;
    
    float Dimension; // Dimension
    float Tinit; // Initial temperature of the annealing simulation
} ;

////// Declarations

// Function that calculates values using the exact D-dimensional random number generator.
std::list<double>* Exact_Generator(SimData &data) ;

// Function that calculates values using the Tsallis random number generator.
std::list<double>* Tsallis_Approximation(SimData &data) ;

// Function that calculates values using the Kleber random number generator.
std::list<double>* Dalligna_Approximation(SimData &data) ;

////// Functions

int main(int argc, char **argv) {
    
    
    //// Variable declaration
    
    std::ifstream InputFile;
    std::ofstream OutputFile;
    
    std::string OrigLine;
    
    // Creates an object that converts between c++ formats, like string<=>int and string<=>float .
    std::stringstream Converter(std::stringstream::in | std::stringstream::out) ;
    
    SimData SimulationData;       // Container of all data necessary for the simulation.
    
    std::list<double> *EG_List = 0, *TA_List = 0, *DA_List = 0; // List of values returned by each number generator.
    
    // Input Processing
    
    try
    {    
        
        InputFile.open( "Input_Info.txt", std::ifstream::in ) ;
        
        while( ! InputFile.eof() )
        {
            std::getline(InputFile, OrigLine) ;
            
            Converter.clear();
            Converter << OrigLine ;
            Converter >> SimulationData.NumLoops ;
            
            std::getline(InputFile, OrigLine) ;
            
            Converter.clear();
            Converter << OrigLine ;
            Converter >> SimulationData.Tinit ;
            
            std::getline(InputFile, OrigLine) ;
            
            Converter.clear();
            Converter << OrigLine ;
            Converter >> SimulationData.qa ;
            
            std::getline(InputFile, OrigLine) ;
            
            Converter.clear();
            Converter << OrigLine ;
            Converter >> SimulationData.qv ;
            
            std::getline(InputFile, OrigLine) ;
            
            Converter.clear();
            Converter << OrigLine ;
            Converter >> SimulationData.qt ;
            
            std::getline(InputFile, OrigLine) ;
            
            Converter.clear();
            Converter << OrigLine ;
            Converter >> SimulationData.myseed ;
        }
        
        // Main processing
        
        std::cout << "\nHello, density!" << std::endl;
        
        std::cout << "GSA parameters" << std::endl;
        std::cout << "Tinit: " << SimulationData.Tinit << "; qa: " << SimulationData.qa << 
        "; qv: " << SimulationData.qv << "; qt: " << SimulationData.qt << "; seed: " << SimulationData.myseed << std::endl;
        
        std::cout << "\nRuning " << SimulationData.NumLoops << " loops for each probability density function... it may take a while." << std::endl; 
        
        
        EG_List = Exact_Generator( SimulationData ) ;
        TA_List = Tsallis_Approximation( SimulationData ) ;
        DA_List = Dalligna_Approximation( SimulationData ) ;
        
        
        // Output
        
        OutputFile.open( "Output_Info.txt", std::ofstream::out ) ;
        
        OutputFile << "GSA parameters" << std::endl;
        OutputFile << "Tinit: " << SimulationData.Tinit << "; qa: " << SimulationData.qa << 
        "; qv: " << SimulationData.qv << "; qt: " << SimulationData.qt << "; seed: " << SimulationData.myseed << "; Total loops: " << SimulationData.NumLoops << std::endl;
        
    //     OutputFile.width(10) ;
    //     OutputFile.right ;
    //     OutputFile.setf(std::ios::fixed,std::ios::floatfield) ;
        OutputFile.precision(8) ;
        
        for ( int loop = 1; loop <= SimulationData.NumLoops; loop++)
        {
            OutputFile << loop << " ; " << EG_List->front() << " ; " << TA_List->front() << " ; " << DA_List->front() << std::endl;
            EG_List->pop_front();
            TA_List->pop_front();
            DA_List->pop_front();
        }
        
    }
    catch ( std::exception &e )
    {
        std::cout << "An error occured (and its probably your fault): " << e.what() << std::endl ;
        
        if ( EG_List != 0 ) delete EG_List;
        
        if ( TA_List != 0 ) delete TA_List;
        
        if ( DA_List != 0 ) delete DA_List;
        
        return 1;
    }
    
    // Clean-up
    
    if ( EG_List != 0 ) delete EG_List;
    
    if ( TA_List != 0 ) delete TA_List;
    
    if ( DA_List != 0 ) delete DA_List;
    
    // THE END
    
    return 0;
}

std::list< double >* Exact_Generator(SimData& data)
{
    //////////////////
    // Exact d-dimensional tsallis random number generator
    // Computer Physics Communications 175 (2006) 708–712
    //////////////////
    
    // Storage of generated values.
    std::list<double> * ReturnList = new std::list<double> ;
    
    /////////////////////////////////////
    // Common variables for the GSA simulation that need initialization
    float qa = data.qa, qt = data.qt, qv = data.qv;
    
    int myseed = data.myseed;
    
    float Dimension = data.Dimension;
    float Tinit = data.Tinit;
    
    // Common variables for the GSA simulation that don't need initialization
    float Tcurr; // Current temperature of the annealing simulation
    
    double Gv; // Visitation distribution value.
    long double Gv_NotNormalized; // Visitation distribution value before normalization.
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Variables for the Temperature calculation.
    long double t_Tinit_times_pow, qt1 ;
    /////////////////////////////////////

    /////////////////////////////////////
    // Variables for the Visitation calculation.
    long double v_p_var ;
    
    long double v_s_var  ;
    long double v_s_sqrt, v_s_13qv ;
    
    long double v_u_var ;
    
    long double v_y_var ;
    
    long double v_LongVar360 ;
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Variables that generate the random numbers.
    int NormMean = 0, NormStdev = 1;
    
    boost::mt19937 Generator ; // Random number generator engine.
    
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Prepare for temperature calculation
    t_Tinit_times_pow = Tinit*(pow(2,(qt - 1)) - 1) ;
    qt1 = qt - 1;
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Variables for the Visitation calculation.
    v_p_var = (3 - qv)/(2*(qv - 1)) ;
    
    v_s_sqrt = sqrt(2 * (qv - 1)) ;
    v_s_13qv = 1/(3 - qv) ;
    
    v_LongVar360 = 360;
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Prepare for random number generation

    Generator.seed( myseed ) ; // Random number generator engine seed.

    boost::normal_distribution<double> NormDist(NormMean,NormStdev);
    
    boost::gamma_distribution<double> GammaDist(v_p_var) ;
    
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > RandNormGenerator(Generator, NormDist) ;
    
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<double> > RandGammaGenerator(Generator, GammaDist) ;
    
    /////////////////////////////////////
    
    
    // Execution loop
    
    std::cout << "Running Exact Generator" << std::endl;
    
    for (int loop = 1; loop <= data.NumLoops; loop++)
    {
        
        Tcurr = t_Tinit_times_pow/( pow((1 + loop),qt1) - 1 ) ; // Temperature T(t)
        
        v_s_var = v_s_sqrt/pow(Tcurr, v_s_13qv) ;
        v_y_var = v_s_var*sqrt( RandGammaGenerator() ) ;
        
        Gv_NotNormalized = RandNormGenerator()/v_y_var;
        
        Gv = fmod( Gv_NotNormalized, v_LongVar360) ;
        
        ReturnList->push_back( Gv );
    }
    
    return ReturnList;
}


std::list< double >* Tsallis_Approximation(SimData& data)
{
    
    //////////////////
    // Function implemented in: Tsallis, Stariolo; Generalized Simulated Annealing (1997)
    //////////////////
    
    // Storage of generated values.
    std::list<double> * ReturnList = new std::list<double> ;
    
    /////////////////////////////////////
    // Common variables for the GSA simulation that need initialization
    float qa = data.qa, qt = data.qt, qv = data.qv;
    
    int myseed = data.myseed;
    
    float Dimension = data.Dimension;
    float Tinit = data.Tinit;
    
    // Common variables for the GSA simulation that don't need initialization
    long double Fact1, Fact2, Fact3, Fact4, Fact5, Fact6, sigmax, nom, den ;
    long double qv1, qv3, sqrtPI, Fact3qv3, qv1qv3 ;
    long double Var360;
    float Tcurr; // Current temperature of the annealing simulation
    
    long double Gv; // Visitation distribution value.
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Variables for the Temperature calculation.
    long double t_Tinit_times_pow, qt1 ;
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Variables that generate the random numbers.
    int NormMean = 0, NormStdev = 1;
    
    boost::mt19937 Generator ; // Random number generator engine.
    
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Prepare for temperature calculation
    t_Tinit_times_pow = Tinit*(pow(2,(qt - 1)) - 1) ;
    qt1 = qt - 1;
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Prepare for random number generation

    Generator.seed( myseed ) ; // Random number generator engine seed.

    boost::normal_distribution<double> NormDist(NormMean,NormStdev);
    
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > RandNormGenerator(Generator, NormDist) ;
    
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Variables for the Visitation calculation.
    Var360 = 360 ;
    
    qv1 = qv - 1 ;
    qv3 = 3.0 - qv ;
    sqrtPI = sqrt(M_PI) ;
    
    qv1qv3 = qv1/qv3 ;
    
    Fact2 = exp((4 - qv)*log(qv1)) ; 
    Fact3 = exp((2 - qv)*log(2)/(qv1)) ;
    
    Fact3qv3 = Fact3*(qv3) ;
    
    Fact5 = ((1.0/(qv1)) - 0.5) ;
    Fact6 = M_PI*(1.0 - Fact5)/sin(M_PI*(1.0 - Fact5))/exp( boost::math::lgamma(2.0 - Fact5) ) ;
    /////////////////////////////////////
    
    // Execution loop
    
    std::cout << "Running Tsallis Approximation" << std::endl;
    
    for (int loop = 1; loop <= data.NumLoops; loop++)
    {
        Tcurr = t_Tinit_times_pow/( pow((1 + loop),qt1) - 1 ) ; // Temperature T(t)
        
        Fact1 = exp(log(Tcurr)/(qv1)) ;
        Fact4 = (sqrtPI*Fact1*Fact2)/Fact3qv3 ;
        sigmax = exp(-qv1qv3*log(Fact6/Fact4)) ;
        
        nom = sigmax*RandNormGenerator() ;
        den = exp(qv1qv3*log(fabs( RandNormGenerator() ))) ;
        
        Gv = fmod( nom/den,Var360) ;
        
        ReturnList->push_back( Gv );
        
    }
    
    return ReturnList;
    
}


std::list< double >* Dalligna_Approximation(SimData& data)
{
    //////////////////
    // Function defined in: Dall'Igna Júnior, Silva, Mundim, Dardenne; Performance and parameterization of the algorithm Simplified Generalized Simulated Annealing (2004)
    //////////////////
    
    // Storage of generated values.
    std::list<double> * ReturnList = new std::list<double> ;
    
    /////////////////////////////////////
    // Common variables for the GSA simulation that need initialization
    float qa = data.qa, qt = data.qt, qv = data.qv;
    
    int myseed = data.myseed;
    
    float Dimension = data.Dimension;
    float Tinit = data.Tinit;
    
    // Common variables for the GSA simulation that don't need initialization
    long double Expo1, Expo2, qv1, Ratio ;
    
    long double Fact1, Fact2, Fact12, Expo3 ;
    
    float Tcurr; // Current temperature of the annealing simulation
    long double Var360;
    long double Gv; // Visitation distribution value.
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Variables for the Temperature calculation.
    long double t_Tinit_times_pow, qt1 ;
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Variables that generate the random numbers.
    
    boost::mt19937 Generator ; // Random number generator engine.
    
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Prepare for temperature calculation
    t_Tinit_times_pow = Tinit*(pow(2,(qt - 1)) - 1) ;
    qt1 = qt - 1;
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Variables for the Visitation calculation.
    Var360 = 360.0;
    
    qv1 = qv - 1.0 ;
    // SGSA
//     Expo1 = (1.0/qv1) - 0.5 ;
    // GSA
    Expo1 = (1.0/qv1) ;
    
    Expo2 = 2.0/(qv - 3.0) ;
    
    // In this case, D = 1
    Fact1 = pow(qv1/M_PI,0.5) ;
    Fact2 = boost::math::tgamma<double>( 1.0/qv1 )/boost::math::tgamma<double>( (1.0/qv1) - 0.5 ) ;
    Fact12 = Fact1*Fact2 ;
    Expo3 = 1.0/(qv - 3.0) ;
    
    /////////////////////////////////////
    
    /////////////////////////////////////
    // Prepare for random number generation

    Generator.seed( myseed ) ; // Random number generator engine seed.

    boost::uniform_01<float> UnifDist;
    
    boost::variate_generator<boost::mt19937&, boost::uniform_01<float> > UnifGenerator(Generator, UnifDist);
    
    /////////////////////////////////////
    
    // Execution loop
    
    std::cout << "Running Dalligna Approximation" << std::endl;
    
    for (int loop = 1; loop <= data.NumLoops; loop++)
    {
        Tcurr = t_Tinit_times_pow/( pow((1 + loop),qt1) - 1 ) ; // Temperature T(t)
        
        //GSA - D=1
        Ratio = qv1*pow(UnifGenerator(),2)/pow(Tcurr,Expo2) ;
        Gv = Fact12*pow(Tcurr,Expo3)/pow(1 + Ratio,Expo1) ;
        
        // SGSA - MUDAR O Expo1
//         Ratio = qv1*pow(UnifGenerator(),2)/pow(Tcurr,Expo2) ;
//         Gv = 1/pow(1 + Ratio,Expo1) ;
        
        // Normalization
        Gv = fmod( Gv,Var360) ;
        
        ReturnList->push_back( Gv );
        
    }
    
    return ReturnList;
    
}
