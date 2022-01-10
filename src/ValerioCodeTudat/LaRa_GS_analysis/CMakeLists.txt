TUDAT_ADD_EXECUTABLE(application_GS_analysis
        "LaRa_GS_analysis.cpp"
        ${Tudat_PROPAGATION_LIBRARIES}
		${Tudat_ESTIMATION_LIBRARIES}
		${Boost_LIBRARIES}
        )