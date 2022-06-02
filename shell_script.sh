# !/bin/bash - screen+bash - run code with source > output_shell.dat
echo "Activate tudat-bundle env"
conda activate tudat-bundle
# Reference RISE+LaRa+PRIDE+corr=0

# # Case 1 (only RISE)
# echo "---CASE 1 (only RISE)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "LaRa_boolean = True -> False"
# sed -i '42s+LaRa_boolean = True+LaRa_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# # Case 2 (only LaRa)
# echo "---CASE 2 (only LaRa PRIDE-False remove-False - corr=0)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "RISE_boolean = False -> True"
# sed -i '41s+RISE_boolean = True+RISE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "PRIDE_boolean = True -> False"
# sed -i '46s+PRIDE_boolean = True+PRIDE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# for var_corr in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99
# do 
	# # Case 3 (only LaRa)
	# echo "---CASE 3 (only LaRa PRIDE-True remove-False - corr=$var_corr)---"
	# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
	# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "Change build path"
	# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "RISE_boolean = False -> True"
	# sed -i '41s+RISE_boolean = True+RISE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "correlation = 0 -> $var_corr"
	# sed -i "52s+correlation = 0+correlation = $var_corr+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "Run code"
	# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo ""
# done

# # Case 4 (RISE + LaRa)
# echo "---CASE 4 (RISE+LaRa PRIDE-False remove-False - corr=0)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "PRIDE_boolean = True -> False"
# sed -i '46s+PRIDE_boolean = True+PRIDE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# for var_corr in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99
# do 
	# # Case 5 (RISE + LaRa)
	# echo "---CASE 5 (RISE+LaRa PRIDE-True remove-False - corr=$var_corr)---"
	# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
	# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "Change build path"
	# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "correlation = 0 -> $var_corr"
	# sed -i "52s+correlation = 0+correlation = $var_corr+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "Run code"
	# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo ""
# done

# lat_value="4.5"

# # Case 6 (RISE + LaRa - latitude changed)
# echo "---CASE 6 (RISE+LaRa PRIDE-False remove-False - corr=0 - latLaRa = $lat_value)---"
# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change build path"
# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change output folder"
# sed -i "68s+POD_RISE+PODlat${lat_value}_RISE+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Change LaRa_reflector_latitude_deg = 18.3 -> $lat_value"
# sed -i "104s+LaRa_reflector_latitude_deg = 18.3+LaRa_reflector_latitude_deg = $lat_value+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "PRIDE_boolean = True -> False"
# sed -i '46s+PRIDE_boolean = True+PRIDE_boolean = False+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo "Run code"
# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
# echo ""

# for var_corr in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 # Add the rest of the var_corr
# do 
	# # Case 7 (RISE + LaRa)
	# echo "---CASE 7 (RISE+LaRa PRIDE-True remove-False - corr=$var_corr - latLaRa = $lat_value)---"
	# echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
	# cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "Change build path"
	# sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "Change output folder"
	# sed -i "68s+POD_RISE+PODlat${lat_value}_RISE+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "Change LaRa_reflector_latitude_deg = 18.3 -> $lat_value"
	# sed -i "104s+LaRa_reflector_latitude_deg = 18.3+LaRa_reflector_latitude_deg = $lat_value+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "correlation = 0 -> $var_corr"
	# sed -i "52s+correlation = 0+correlation = $var_corr+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo "Run code"
	# python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
	# echo ""
# done

# noise_factor="10"

# for var_corr in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 # Add the rest of the var_corr
# do
        # # Case 8 (RISE + LaRa)
        # echo "---CASE 8 (RISE+LaRa PRIDE-True remove-False -corr=$varr_corr - DSN-LaRa noise reduced by a factor of $noise_factor)"
        # echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py"
        # cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
        # echo "Change build path"
        # sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
        # echo "Change output folder"
        # sed -i "68s+POD_RISE+PODnoise_RISE+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
        # echo "correlation = 0 -> $var_corr"
        # sed -i "52s+correlation = 0+correlation = $var_corr+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
        # echo "Change noise for DSN-LaRa"
        # sed -i "702s+)))+)/$noise_factor))+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
        # echo "Run code"
        # python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
        # echo ""
# done	

#for var_corr in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#do 
#	for station_n in 1 2 3 4 5 6 7 8 9 10
#	do
#		# Case 9 (only LaRa) - V&V
#		echo "---CASE 9 (only LaRa PRIDE-True remove-False -corr=$varr_corr - stations $station_n)"
#		echo "Copy py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_GS_ver.py"
#		cp -r ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_GS_ver.py ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
#		echo "Change build path"
#		sed -i '14s+sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")+sys.path.insert(0, "/home/carlos/tudat-bundle/build_release/tudatpy/")+' ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
#		echo "Change correlation = 0 -> $var_corr"
#		sed -i "52s+correlation = 0+correlation = $var_corr+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
#		echo "Change number of stations = 1 -> $station_n"
#		sed -i "67s+receiving_station_number = 1+receiving_station_number = $station_n+" ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
#		echo "Run code"
#		python ./src/py_preciseorbitdetermination_InSight_LaRa_DSN_PRIDE_server.py
#		echo ""
#	done
#done