/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120522    A. Ronse          First creation of code.
 *
 *    References
 *      Williams, Dr. David R., "Moon Fact Sheet", NASA (National Space Science Data Center),
 *         http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html, last accessed: 22 May 2012
 *
 *    Notes
 *
 */

#include <tudat/simulation/estimation.h>
#include <tudat/astro/orbit_determination/estimatable_parameters/periodicSpinVariation.h>

int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::observation_models;
    using namespace tudat::orbit_determination;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::ephemerides;
    using namespace tudat::propagators;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::coordinate_conversions;
    using namespace tudat::ground_stations;
    using namespace tudat::observation_models;
    using namespace tudat::statistics;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;

    bodyNames.push_back( "Saturn" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Venus" );
    bodyNames.push_back( "Mercury" );
    bodyNames.push_back( "Sun" );

    // Specify initial and final time

    double initialEphemerisTime = ( 2459215.5 - JULIAN_DAY_ON_J2000 ) *86400; // 1/01/2021 00:00:00
    int numberOfSimulationDays = 700; // le Maistre simulation
    double finalEphemerisTime = initialEphemerisTime + numberOfSimulationDays * physical_constants::JULIAN_DAY;

    // Create bodies needed in simulation
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings(bodyNames, initialEphemerisTime - physical_constants::JULIAN_DAY, finalEphemerisTime + physical_constants::JULIAN_DAY);

    bodySettings[ "Moon" ]->ephemerisSettings->resetFrameOrigin( "Sun" );

    bodySettings[ "Mars" ]->rotationModelSettings = getHighAccuracyMarsRotationModel( initialEphemerisTime, finalEphemerisTime );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE GROUND STATIONS AND LANDER      ////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create transmitter ground stations from GroundStationDatabase.
    std::vector< std::string > groundStationNames;
//    groundStationNames.push_back( "DSS14" );
    groundStationNames.push_back( "DSS63" );
//    groundStationNames.push_back( "TIDBIN64" ); //DSS43

    // Creat PRIDE ground stations from GroundStationDatabase.
    std::ifstream readFile (".src/gs_location.txt");

    int num_lines = 0;

    if ( readFile.is_open( ) )
    {
        std::string line;

        while ( getline(readFile, line) )
        {
            groundStationNames.push_back( line );

            ++num_lines;
        }

    }
    else
    {
        std::cerr << "Couldn't open gs_location.txt for reading.\n";
    }


    createGroundStations( bodyMap.at( "Earth" ), groundStationNames );

    createGroundStation( bodyMap.at( "Mars" ), "LaRa",
                         ( Eigen::Vector3d( ) << spice_interface::getAverageRadius("Mars"),
                           unit_conversions::convertDegreesToRadians( 18.20 ),
                           unit_conversions::convertDegreesToRadians( 335.45 )).finished( ),
                         spherical_position);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    SelectedAccelerationMap accelerationMap;

    accelerationMap[ "Mars" ][ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Mars" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Mars" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Mars" ][ "Venus" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Mars" ][ "Mercury" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationMap[ "Mars" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );

    // Define list of bodies to propagate
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Mars" );

    // Define central bodies to use in propagation.
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "SSB" );

    // Create acceleration models and propagation settings.
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Get initial state vector as input to integration.
    Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                bodiesToPropagate, centralBodies, bodyMap, initialEphemerisTime );

    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > TraslationalpropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime );

    double initialTimeStep = 1.0;
    double minimumStepSize = initialTimeStep;
    double maximumStepSize = 60.0;
    double relativeErrorTolerance = 1.0E-14;
    double absoluteErrorTolerance = 1.0E-14;

    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< double > >
            ( initialEphemerisTime, initialTimeStep,
              RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78, minimumStepSize,
              maximumStepSize, relativeErrorTolerance, absoluteErrorTolerance);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE LINK ENDS FOR OBSERVATIONS            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create list of link ends in which station is receiver-transmitter and in which lander is reflector (with lander other link end).
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    std::vector< LinkEnds > twoWayLinkEnds;
    for( unsigned int i = 0; i < groundStationNames.size( ) ; i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", "DSS63" );
        linkEnds[ reflector1 ] = std::make_pair( "Mars", "LaRa" );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        twoWayLinkEnds.push_back( linkEnds );

        // Define link ends to be used for 2-way doppler
        linkEndsPerObservable[ two_way_doppler ].push_back( twoWayLinkEnds[ i ] );
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS THAT ARE TO BE ESTIMATED      ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;

    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Mars", systemInitialState, "SSB" ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", ground_station_position, "LaRa" ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", core_factor ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", free_core_nutation_rate ) );
    parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >( "Mars", periodic_spin_variation ) );
    parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >( "Mars", polar_motion_amplitude ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap, accelerationModelMap );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE OBSERVATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create light-time correction settings.
    std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections;
    std::vector< std::string > perturbingBodies;
    perturbingBodies.push_back("Sun");
    lightTimeCorrections = std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                perturbingBodies );

    // Define Uplink OneWay Doppler Observation Settings
    std::shared_ptr< DirectFirstOrderDopplerProperTimeRateSettings > uplinktransmitterProperTimeRateSettings =
            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" );
    std::shared_ptr< DirectFirstOrderDopplerProperTimeRateSettings > uplinkreceiverProperTimeRateSettings =
            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" );
    std::shared_ptr< OneWayDopplerObservationSettings > uplinkOneWayDopplerSettings =
            std::make_shared< OneWayDopplerObservationSettings >(
                lightTimeCorrections, uplinktransmitterProperTimeRateSettings, uplinkreceiverProperTimeRateSettings );

    // Define Downlink OneWay Doppler Observation Settings
    std::shared_ptr< DirectFirstOrderDopplerProperTimeRateSettings > downlinktransmitterProperTimeRateSettings =
            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" );
    std::shared_ptr< DirectFirstOrderDopplerProperTimeRateSettings > downlinkreceiverProperTimeRateSettings =
            std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" );
    std::shared_ptr< OneWayDopplerObservationSettings > downlinkOneWayDopplerSettings =
            std::make_shared< OneWayDopplerObservationSettings >(
                lightTimeCorrections, downlinktransmitterProperTimeRateSettings, downlinkreceiverProperTimeRateSettings );

    // Define TwoWay Doppler Observation Settings
    std::shared_ptr< TwoWayDopplerObservationSettings > observationSettings =
            std::make_shared< TwoWayDopplerObservationSettings >(
                uplinkOneWayDopplerSettings, downlinkOneWayDopplerSettings );

    // Iterate over all observable types and associated link ends, and creating settings for observation
    observation_models::ObservationSettingsMap observationSettingsMap;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        // Define ObservationSettingsMap
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsMap.insert( std::make_pair( currentLinkEndsList.at( i ), observationSettings ) );
        }
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          INITIALIZE ORBIT DETERMINATION OBJECT     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create orbit determination object (propagate orbit, create observation models)
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodyMap, parametersToEstimate, observationSettingsMap,
                integratorSettings, TraslationalpropagatorSettings );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          SIMULATE OBSERVATIONS                     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Define time of first observation
    double observationTimeStart = initialEphemerisTime + physical_constants::JULIAN_DAY;

    // Define time between two observations
    double  observationInterval = 60.0;

    // Define numbers of weeks
    double numberOfSimulationWeeks = numberOfSimulationDays / 7.0;

    // Simulate observations for each week in simulation
    std::vector< double > baseTimeList;
    for( int i = 0; i < numberOfSimulationWeeks ; i++ )
    {
        for( unsigned int j = 0; j < physical_constants::JULIAN_DAY / observationInterval; j++ )
        {
            for( int obsPerWeek = 0; obsPerWeek < 2; obsPerWeek++ )
            {

            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 7.0 * physical_constants::JULIAN_DAY +
                                    static_cast< double >( obsPerWeek ) * 3.25 * physical_constants::JULIAN_DAY +
                                    static_cast< double >( j ) * observationInterval );
            }
        }
    }

    // Create measureement simulation input
    std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< ObservationSimulationTimeSettings< double > > > >
            measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        // Define observable type and link ends
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        // Define observation times and reference link ends
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                    std::make_shared< TabulatedObservationSimulationTimeSettings< double > >( transmitter, baseTimeList );
        }
    }

    // Create observation viability settings and calculators
    std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettings;
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                minimum_elevation_angle, std::make_pair( "Earth", "" ), "",
                                                unit_conversions::convertDegreesToRadians( 20.0 ) ) );
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                minimum_elevation_angle, std::make_pair( "Mars", "" ), "",
                                                unit_conversions::convertDegreesToRadians( 35.0 ) ) );
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                maximum_elevation_angle, std::make_pair( "Mars", "" ), "",
                                                unit_conversions::convertDegreesToRadians( 45.0 ) ) );
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                body_avoidance_angle, std::make_pair( "Earth", "" ), "Sun",
                                                unit_conversions::convertDegreesToRadians( 20.0 ) ) );
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                body_occultation, std::make_pair( "Earth", "" ), "Moon" ) );
    PerObservableObservationViabilityCalculatorList viabilityCalculators = createObservationViabilityCalculators(
                bodyMap, linkEndsPerObservable, observationViabilitySettings );

    // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with
    // reference link ends.
    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
            SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Define noise levels
    double dopplerNoise = 0.05E-3 / physical_constants::SPEED_OF_LIGHT ; // Doppler error budget Denhant et all (2009)

    // Create noise functions per observable
    std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;

    noiseFunctions[ two_way_doppler ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                           normal_boost_distribution, { 0.0, dopplerNoise }, 0.0 ), std::placeholders::_1 );

    // Simulate observations
    PodInputDataType observationsAndTimes = simulateObservationsWithNoise< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), noiseFunctions,
                viabilityCalculators );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////    PERTURB PARAMETER VECTOR AND ESTIMATE PARAMETERS     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );
//    parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 100.0 ); // Mars position [m]
//    parameterPerturbation.segment(  3, 3 ) = Eigen::Vector3d::Constant( 1.0 ); // Mars velocity [m/s]
//    parameterPerturbation( 6 ) =  0 ; // core factor of the celestial body of (Mars) [rad]
//    parameterPerturbation( 7 ) =  0 ; // free core nutation rate of the celestial bodyof (Mars) [rad/s]
    parameterPerturbation( 8 ) =  19E3 ; // Lander x-position [m]
    parameterPerturbation( 9 ) =  159E3 ; // Lander y-position [m]
    parameterPerturbation( 10 ) =  53E3 ; // Lander z-position [m]



    initialParameterEstimate += parameterPerturbation;

    // Define a priori covariance
    Eigen::MatrixXd InverseAPriopriCovariance =
            Eigen::MatrixXd::Zero( initialParameterEstimate.rows(), initialParameterEstimate.rows() );

//    for( unsigned int i = 0; i <= 5; i++ )
//    {
//        InverseAPriopriCovariance( i, i ) = 1.0 ;
//    }
//    InverseAPriopriCovariance( 6, 6 ) =  1.0 / std::pow( 200E3, 2 );
//    InverseAPriopriCovariance( 7, 7 ) =  1.0 / std::pow( 200E3, 2 );
//    InverseAPriopriCovariance( 8, 8 ) =  1.0 / std::pow( 200E3, 2 );

    std::shared_ptr< PodInput< double, double > > podInput = std::make_shared< PodInput< double, double > >(
                observationsAndTimes, ( initialParameterEstimate ).rows( ) );

//    // Define estimation input
//    std::shared_ptr< PodInput< double, double > > podInput =
//            std::make_shared< PodInput< double, double > >(
//                observationsAndTimes, initialParameterEstimate.rows( ),
//                 InverseAPriopriCovariance, initialParameterEstimate - truthParameters);

//    // Define observation weights (constant per observable type)
//    std::map< observation_models::ObservableType, double > weightPerObservable;
//    weightPerObservable[ two_way_doppler ] = 1.0 / ( dopplerNoise * dopplerNoise );
//    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );

    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( 1 ) );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::string outputFolder = ".src/";

    // Print true estimation error, limited mostly by numerical error
    Eigen::VectorXd estimationError = podOutput->parameterEstimate_ - truthParameters;

    std::cout << "True estimation error is:   " << std::endl << ( estimationError ).transpose( ) << std::endl;
    std::cout << "Formal estimation error is: " << std::endl << podOutput->getFormalErrorVector( ).transpose( ) << std::endl;
    std::cout << "True to form estimation error ratio is: " << std::endl <<
                 ( podOutput->getFormalErrorVector( ).cwiseQuotient( estimationError ) ).transpose( ) << std::endl;

    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                     "EstimationInformationMatrix.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_,
                                     "EstimationInformationMatrixNormalization.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
                                     "ObservationTimes.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
                                     "ObservationLinkEnds.dat", 16,
                                     outputFolder );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;

}

