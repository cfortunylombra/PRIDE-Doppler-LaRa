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
#include <tudat/io/basicInputOutput.h>
#include <iostream>
#include <fstream>

int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    using namespace tudat;
    using namespace tudat::interpolators;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::ephemerides;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::coordinate_conversions;
    using namespace tudat::ground_stations;
    using namespace tudat::unit_conversions;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////     CREATE ENVIRONMENT       /////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Load spice kernels.
    loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );

    // Specify initial and final time
//    double initialEphemerisTime = ( 2459252.5 - JULIAN_DAY_ON_J2000 ) *86400; // 7/02/2021 00:00:00
    //double initialEphemerisTime = ( 2459274.5 - JULIAN_DAY_ON_J2000 ) *86400; // 1/03/2021 00:00:00
    double initialEphemerisTime = ( 2460004.5 - JULIAN_DAY_ON_J2000 ) *86400; // 01/03/2023 00:00:00
//    int numberOfSimulationDays = 365; // One Earth year
//    int numberOfSimulationDays = 687; // One Martian year
    int numberOfSimulationDays = 700; // le Maistre simulation
    double finalEphemerisTime = initialEphemerisTime + numberOfSimulationDays * physical_constants::JULIAN_DAY;

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames,
                                                            initialEphemerisTime-physical_constants::JULIAN_DAY,
                                                            finalEphemerisTime+physical_constants::JULIAN_DAY,
                                                            "SSB", "ECLIPJ2000", 60 );

    bodySettings.get( "Mars" )->rotationModelSettings = simulation_setup::getHighAccuracyMarsRotationModel();

    SystemOfBodies bodyMap = createSystemOfBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE GROUND STATIONS AND LANDER      ////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;
//    groundStationNames.push_back( "DSS14" );
    groundStationNames.push_back( "DSS63" );

//    std::ifstream readFile ("/home/cfortunylombra/tudat-bundle/tudat/examples/tudat/LaRa_GS_analysis/src/gs_location.txt");
//    std::ifstream readFile ("/Users/valeriofilice/Tudat/tudatBundle/tudatApplications/Outputs/LaRa_GS_analysis_EVN/EVN_location.txt");
//    std::ifstream readFile ("/Users/valeriofilice/Tudat/tudatBundle/tudatApplications/Outputs/LaRa_GS_analysis_VLA/VLA_location.txt");
//    std::ifstream readFile ("/Users/valeriofilice/Tudat/tudatBundle/tudatApplications/Outputs/LaRa_GS_analysis_Extra/Extra_location.txt");
//    std::ifstream readFile ("/Users/valeriofilice/Tudat/tudatBundle/tudatApplications/Outputs/LaRa_GS_analysis_ATNF/ATNF_location.txt");
//    std::ifstream readFile ("/Users/valeriofilice/Tudat/tudatBundle/tudatApplications/Outputs/LaRa_GS_analysis_DSN/DSN_station.txt.txt");
//    int num_lines = 0;
//
//    if ( readFile.is_open( ) )
//    {
//        std::string line;
//
//        while ( getline(readFile, line) )
//        {
//            groundStationNames.push_back( line );
//
//            ++num_lines;
//        }
//
//    }
//    else
//    {
//        std::cerr << "Couldn't open locations.dat.txt for reading.\n";
//    }

//    std::cout << num_lines << std::endl;

    //createGroundStations( bodyMap.at( "Earth" ), groundStationNames );

//    createGroundStation( bodyMap.at( "Mars" ), "LaRa",
//                         ( Eigen::Vector3d( ) << spice_interface::getAverageRadius("Mars"),
//                           unit_conversions::convertDegreesToRadians( -2.0 ),
//                           unit_conversions::convertDegreesToRadians( 354.0 )).finished( ),
//                         spherical_position);

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
                           unit_conversions::convertDegreesToRadians( 18.3 ),
                           unit_conversions::convertDegreesToRadians( 335.37 )).finished( ),
                         spherical_position);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////             GROUND STATIONS ELEVATION HISTORY            //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define time of first observation
    double observationTimeStart = initialEphemerisTime + physical_constants::JULIAN_DAY;

    // Define time between two observations
    double  observationInterval = 60.0;

    // Define numbers of weeks
    double numberOfSimulationWeeks = numberOfSimulationDays / 7.0;

    // Simulate observations for each week in simulation
    std::vector< double > baseTimeList;
    for( int i = 0; i < numberOfSimulationWeeks; i++ )
    {
        for( unsigned int j = 0; j < physical_constants::JULIAN_DAY / observationInterval; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 7.0 * physical_constants::JULIAN_DAY +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    std::shared_ptr< GroundStation > LaRa = bodyMap.at( "Mars" )->getGroundStation( "LaRa" );
    std::shared_ptr< GroundStationState > LaRaNominalStateObject = LaRa->getNominalStationState( );
    std::shared_ptr< PointingAnglesCalculator > LaRaPointingAngleCalculatorObject =
            LaRa->getPointingAnglesCalculator( );

    std::shared_ptr< GroundStation > DSS63 = bodyMap.at( "Earth" )->getGroundStation( "DSS63" );
    std::shared_ptr< GroundStationState > DSS63NominalStateObject = DSS63->getNominalStationState( );
    std::shared_ptr< PointingAnglesCalculator > DSS63PointingAngleCalculatorObject =
            DSS63->getPointingAnglesCalculator( );

    std::vector< double > observationTimes;
    std::vector< double > earthElevation;
    std::vector< double > earthAzimuth;
    std::vector< double > DSS63ObservationTime;
    std::vector < double > DSS63Elevation;

    for ( double currentTime : baseTimeList )
    {
        observationTimes.push_back( currentTime );
        Eigen::Matrix3d rotationFromMarsLocalFrametoInertialFrame =bodyMap.at( "Mars" )->getRotationalEphemeris( )->getRotationMatrixToBaseFrame(currentTime);
        Eigen::Matrix3d rotationFromEarthLocalFrametoInertialFrame = bodyMap.at( "Earth" )->getRotationalEphemeris( )->getRotationMatrixToBaseFrame(currentTime);

//        earthElevation.push_back( LaRaPointingAngleCalculatorObject->calculateElevationAngle(
//                                  ( bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( currentTime ) -
//                                  bodyMap[ "Mars" ]->getStateInBaseFrameFromEphemeris( currentTime ) ).segment( 0, 3 ), currentTime ) );
//        earthAzimuth.push_back( LaRaPointingAngleCalculatorObject->calculationAzimuthAngle(
//                                    ( bodyMap[ "Earth" ]->getStateInBaseFrameFromEphemeris( currentTime ) -
//                                    bodyMap[ "Mars" ]->getStateInBaseFrameFromEphemeris( currentTime ) ).segment( 0, 3 ), currentTime ) );

        earthElevation.push_back(
                    LaRaPointingAngleCalculatorObject->calculateElevationAngle(
                        ( bodyMap.at("Earth")->getStateInBaseFrameFromEphemeris( currentTime ) ).segment( 0, 3 ) -
                        rotationFromMarsLocalFrametoInertialFrame * LaRaNominalStateObject->getCartesianPositionInTime( currentTime ),
                        currentTime ) );
        earthAzimuth.push_back(
                    LaRaPointingAngleCalculatorObject->calculateAzimuthAngle(
                        ( bodyMap.at("Earth")->getStateInBaseFrameFromEphemeris( currentTime ) ).segment( 0, 3 ) -
                        rotationFromMarsLocalFrametoInertialFrame * LaRaNominalStateObject->getCartesianPositionInTime( currentTime ) ,
                        currentTime ) );

        if ( ( earthElevation.back( ) >= unit_conversions::convertDegreesToRadians( 35.0 ) ) &&
             ( earthElevation.back( ) <= unit_conversions::convertDegreesToRadians( 45.0 ) ) ) // LaRa antenna constrain
        {
            DSS63ObservationTime.push_back( currentTime );

            DSS63Elevation.push_back(
                        DSS63PointingAngleCalculatorObject->calculateElevationAngle(
//                            rotationFromMarsLocalFrametoInertialFrame * LaRaNominalStateObject->getCartesianPositionInTime( currentTime ) -
                                bodyMap.at( "Mars" )->getStateInBaseFrameFromEphemeris( currentTime ).segment( 0, 3 ) -
                            rotationFromEarthLocalFrametoInertialFrame * DSS63NominalStateObject->getCartesianPositionInTime( currentTime ),
                            currentTime ));

        }

    }

    std::map< std::string, std::shared_ptr< GroundStation > > GroundStationMap =
            bodyMap.at( "Earth" )->getGroundStationMap( );

    std::vector < int > groundStationIDs;
    std::vector< double > groundStationObservationTimes;
    std::vector< double > groundStationElevations;
    std::vector< string > groundStationName;
    int Id = 0;

    for( std::map< std::string, std::shared_ptr< GroundStation > >::iterator gsIterator = GroundStationMap.begin( );
        gsIterator != GroundStationMap.end( ); gsIterator++ )
    {
        std::shared_ptr< GroundStation > currentGroudStation = gsIterator->second;

        if ( currentGroudStation->getStationId() == "DSS63" ) continue;

        std::shared_ptr< GroundStationState > currentNominalStateObject = currentGroudStation->getNominalStationState( );
        std::shared_ptr< PointingAnglesCalculator > groudStationPointingAngleCalculatorObject =
                currentGroudStation->getPointingAnglesCalculator( );

        int ind = 0;
        for ( double currentTime : DSS63ObservationTime )
        {
            if ( DSS63Elevation[ ind ] >= unit_conversions::convertDegreesToRadians( 20.0 ) )
            {
                groundStationObservationTimes.push_back( currentTime );
                Eigen::Matrix3d rotationFromEarthLocalFrametoInertialFrame = bodyMap.at( "Earth" )->getRotationalEphemeris( )->getRotationMatrixToBaseFrame(currentTime);

                groundStationElevations.push_back( groudStationPointingAngleCalculatorObject->calculateElevationAngle(
//                               rotationFromMarsLocalFrametoInertialFrame *
//                                                       LaRaNominalStateObject->getCartesianPositionInTime( currentTime ) -
                        bodyMap.at( "Mars" )->getStateInBaseFrameFromEphemeris( currentTime ).segment( 0, 3 ) -
                       rotationFromEarthLocalFrametoInertialFrame * currentNominalStateObject->getCartesianPositionInTime( currentTime ) ,
                       currentTime ) );

                groundStationIDs.push_back( Id );

            }

            ++ind;
        }

        Id += 1;

    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::string outputFolder = "/home/cfortunylombra/tudat-bundle/tudat/examples/tudat/LaRa_GS_analysis/output/";

//    std::string outputSubFolder = "LaRa_GS_analysis_EVN";
//    std::string outputSubFolder = "LaRa_GS_analysis_VLA/";
//    std::string outputSubFolder = "LaRa_GS_analysis_Extra/";
//    std::string outputSubFolder = "LaRa_GS_analysis_DSN/";
//    std::string outputSubFolder = "LaRa_GS_analysis_ATNF/";

    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( observationTimes ),
                                     "time.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( DSS63ObservationTime ),
                                     "DSS63ObservationTime.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( DSS63Elevation ),
                                     "DSS63Elevation.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( earthElevation ),
                                     "earthElevations.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( earthAzimuth ),
                                     "earthAzimuth.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( groundStationObservationTimes ),
                                     "groundStationObservationTimes.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( groundStationElevations ),
                                     "groundStationElevations.dat", 16,
                                     outputFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( groundStationIDs ),
                                     "groundStationIDs.dat", 16,
                                     outputFolder );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
