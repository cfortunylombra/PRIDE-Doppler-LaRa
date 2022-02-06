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

    std::cout.precision(15);
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
    double initialEphemerisTime = ( 2459215.5 - JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;  // 1/01/2021 00:00:00
    double numberOfSimulationDays = 49.0; // le Maistre simulation
    double finalEphemerisTime = initialEphemerisTime + numberOfSimulationDays * physical_constants::JULIAN_DAY;

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames,
                                                            initialEphemerisTime-physical_constants::JULIAN_DAY,
                                                            finalEphemerisTime+physical_constants::JULIAN_DAY,
                                                            "SSB", "ECLIPJ2000", 60 );

    bodySettings.get( "Moon" )->ephemerisSettings->resetFrameOrigin( "Sun" );

    bodySettings.get( "Mars" )->rotationModelSettings = simulation_setup::getHighAccuracyMarsRotationModel();

    SystemOfBodies bodyMap = createSystemOfBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE GROUND STATIONS AND LANDER      ////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create transmitter ground stations from GroundStationDatabase.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "BADARY" );
    groundStationNames.push_back( "CEDUNA" );
    groundStationNames.push_back( "HARTRAO" );
    groundStationNames.push_back( "HART15M" );
    groundStationNames.push_back( "HOBART12" );
    groundStationNames.push_back( "HOBART26" );
    groundStationNames.push_back( "TIANMA65" );
    groundStationNames.push_back( "WARK30M" );
    groundStationNames.push_back( "EFLSBERG" );
    groundStationNames.push_back( "IRBENE" );
    groundStationNames.push_back( "YEBES40M" );
    groundStationNames.push_back( "MEDICINA" );
    groundStationNames.push_back( "WETTZELL" );
    groundStationNames.push_back( "ONSALA60" );
    groundStationNames.push_back( "WRT0" );

    createGroundStation(bodyMap.at("Earth"), "DSS63",
                        (Eigen::Vector3d() << 4849092.6814,-360180.5350,4115109.1298). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "BADARY",
                        (Eigen::Vector3d() <<-838201.2618,3865751.5589,4987670.8708). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "CEDUNA",
                        (Eigen::Vector3d() <<-3753442.7457,3912709.7530,-3348067.6095). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "HARTRAO",
                        (Eigen::Vector3d() <<5085442.7721,2668263.9300,-2768696.6299). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "HART15M",
                        (Eigen::Vector3d() <<5085490.8071,2668161.6274,-2768692.5007). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "HOBART12",
                        (Eigen::Vector3d() <<-3949991.0556,2522421.2681,-4311707.7596). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "HOBART26",
                        (Eigen::Vector3d() <<-3950237.6192,2522347.7349,-4311561.5974). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "TIANMA65",
                        (Eigen::Vector3d() <<-2826708.8081,4679236.9722,3274667.4495). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "WARK30M",
                        (Eigen::Vector3d() <<-5115423.680,477880.102,-3767040.597). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "EFLSBERG",
                        (Eigen::Vector3d() <<4033947.1525,486990.8961,4900431.0604). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "IRBENE",
                        (Eigen::Vector3d() <<3183649.341,1276902.985,5359264.715). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "YEBES40M",
                        (Eigen::Vector3d() <<4848761.7579,-261484.0570,4123085.1343). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "MEDICINA",
                        (Eigen::Vector3d() <<4461369.5682,919597.2489,4449559.4702). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "WETTZELL",
                        (Eigen::Vector3d() <<4075539.5173,931735.6497,4801629.6028). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "ONSALA60",
                        (Eigen::Vector3d() <<3370605.7035,711917.8146,5349830.9852). finished( ),
                        cartesian_position);

    createGroundStation(bodyMap.at("Earth"), "WRT0",
                        (Eigen::Vector3d() <<3828767.1338,442446.1588,5064921.5700). finished( ),
                        cartesian_position);


    createGroundStation( bodyMap.at( "Mars" ), "LaRa",
                         ( Eigen::Vector3d( ) << spice_interface::getAverageRadius("Mars"),
                           unit_conversions::convertDegreesToRadians( 18.20 ),
                           unit_conversions::convertDegreesToRadians( 335.45 )).finished( ),
                         spherical_position);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    SelectedAccelerationMap accelerationMap;

    accelerationMap[ "Mars" ][ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Mars" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Mars" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Mars" ][ "Venus" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Mars" ][ "Mercury" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Mars" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

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

    // Define propagation termination conditions
    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( finalEphemerisTime );

    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings );

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
            createParametersToEstimate< double >( parameterNames, bodyMap, propagatorSettings );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE OBSERVATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Iterate over all observable types and associated link ends, and creating settings for observation
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;

    // Define ObservationSettingsMap
    for( unsigned int i = 0; i < twoWayLinkEnds.size( ); i++ )
    {
        std::shared_ptr< TwoWayDopplerObservationSettings > twoWayObservationSettings =
                std::make_shared< TwoWayDopplerObservationSettings >(
                        twoWayLinkEnds[i] );

        observationSettingsList.push_back( twoWayObservationSettings );

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          INITIALIZE ORBIT DETERMINATION OBJECT     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create orbit determination object (propagate orbit, create observation models)
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodyMap, parametersToEstimate, observationSettingsList,
                integratorSettings, propagatorSettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          SIMULATE OBSERVATIONS                     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Define time of first observation
    double observationTimeStart = initialEphemerisTime + physical_constants::JULIAN_DAY;

    // Define time between two observations
    double  observationInterval = 60.0*60.0;

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
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        // Define observable type and link ends
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        // Define observation times and reference link ends
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            currentObservable, currentLinkEndsList.at( i ), baseTimeList, transmitter ) );
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
    //observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
    //                                            maximum_elevation_angle, std::make_pair( "Mars", "" ), "",
    //                                            unit_conversions::convertDegreesToRadians( 45.0 ) ) );
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                body_avoidance_angle, std::make_pair( "Earth", "" ), "Sun",
                                                unit_conversions::convertDegreesToRadians( 20.0 ) ) );
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                body_occultation, std::make_pair( "Earth", "" ), "Moon" ) );

    //addViabilityToObservationSimulationSettings(
    //        measurementSimulationInput,
    //        observationViabilitySettings );


    // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with
    // reference link ends.
    //typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    //typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
    //        SingleObservablePodInputType;
    //typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Define noise levels
    double dopplerNoise = 0.05E-3 / physical_constants::SPEED_OF_LIGHT ; // Doppler error budget Denhant et all (2009)

    // Create noise functions per observable
    //std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;

    //noiseFunctions[ two_way_doppler ] =
    //        std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
    //                   createBoostContinuousRandomVariableGeneratorFunction(
    //                       normal_boost_distribution, { 0.0, dopplerNoise }, 0.0 ), std::placeholders::_1 );

    // Simulate observations
    std::shared_ptr< observation_models::ObservationCollection< double, double > > observationsAndTimes = simulateObservations< >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodyMap );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////    PERTURB PARAMETER VECTOR AND ESTIMATE PARAMETERS     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );

    std::cout << "Initial parameter estimate is: " << std::endl << (initialParameterEstimate).transpose() << std::endl;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );
    parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1000.0 ); // Mars position [m]
    parameterPerturbation.segment(  3, 3 ) = Eigen::Vector3d::Constant( 10.0 ); // Mars velocity [m/s]
    parameterPerturbation( 6 ) =  1E-2 ;
    parameterPerturbation( 7 ) =  1E-7 ;
    parameterPerturbation.segment(  8, 3 ) = Eigen::Vector3d::Constant( 100.0 );
    parameterPerturbation.segment(  11, 28 ) = Eigen::Vector3d::Constant( 1E-9 );

    std::cout << "Perturbation vector is: " << std::endl << (parameterPerturbation).transpose() << std::endl;

    // Define a priori covariance
    //Eigen::MatrixXd InverseAPriopriCovariance =
    //       Eigen::MatrixXd::Zero( initialParameterEstimate.rows(), initialParameterEstimate.rows() );

    //    for( unsigned int i = 0; i <= 5; i++ )
    //    {
    //        InverseAPriopriCovariance( i, i ) = 1.0 ;
    //    }
    //    InverseAPriopriCovariance( 6, 6 ) =  1.0 / std::pow( 200E3, 2 );
    //    InverseAPriopriCovariance( 7, 7 ) =  1.0 / std::pow( 200E3, 2 );
    //    InverseAPriopriCovariance( 8, 8 ) =  1.0 / std::pow( 200E3, 2 );

    std::shared_ptr< PodInput< double, double > > podInput = std::make_shared< PodInput< double, double > >(
                observationsAndTimes, ( initialParameterEstimate ).rows( ), Eigen::MatrixXd::Zero( 0, 0 ), parameterPerturbation);

    //    // Define observation weights (constant per observable type)
    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ two_way_doppler ] = 1.0 / ( dopplerNoise * dopplerNoise );
    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );

    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput);//, std::make_shared< EstimationConvergenceChecker >( 1 ) );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::string outputFolder = "/home/cfortunylombra/tudat-bundle/tudat/examples/tudat/LaRa_POD/src/";

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
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( observationsAndTimes->getConcatenatedTimeVector( ) ),
                                     "ObservationTimes.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(observationsAndTimes->getConcatenatedLinkEndIds()),
                                     "ObservationLinkEnds.dat", 16,
                                     outputFolder );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;

}
