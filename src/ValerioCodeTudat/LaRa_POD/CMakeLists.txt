TUDAT_ADD_EXECUTABLE(application_LaRa
        "LaRa_POD.cpp"
        ${Tudat_PROPAGATION_LIBRARIES}
	${Tudat_ESTIMATION_LIBRARIES}
	${Boost_LIBRARIES}
        )