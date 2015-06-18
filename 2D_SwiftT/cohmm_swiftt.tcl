namespace eval cohmm_swiftt {
	namespace export cohmm_swiftt
	
	set version 1.0
	set MyDescription "2D CoHMM DaD"

	proc initWrapper {doKriging doCoMD bDims bDt bDelta bGamma} {
		#Unpack blobs into lists
		set deltaPtr [ lindex $bDelta 0]
		set dimsPtr [ lindex $bDims 0 ]
		set dtPtr [ lindex $bDt 0 ]
		set gammaPtr [ lindex $bGamma 0]
		#Convert ptrs to swig ptrs
		set dtPtr [ blobutils_cast_int_to_dbl_ptr $dtPtr]
		set deltaPtr [ blobutils_cast_int_to_dbl_ptr $deltaPtr]
		set gammaPtr [ blobutils_cast_int_to_dbl_ptr $gammaPtr]
		set dimsPtr [ blobutils_cast_int_to_int_ptr $dimsPtr]
		#Load library
		load ../bld/libcohmm_swiftt.so
		#Call function
		set retValue [ initEverything $doKriging $doCoMD $dimsPtr $dtPtr $deltaPtr $gammaPtr ]
		#Return retValue
		return $retValue
	}

	proc prepFirstWrapper {doKriging doCoMD bDims bDt bDelta bGamma step} {
		#Unpack blobs into lists
		set deltaPtr [ lindex $bDelta 0]
		set dimsPtr [ lindex $bDims 0 ]
		set dtPtr [ lindex $bDt 0 ]
		set gammaPtr [ lindex $bGamma 0]
		#Convert ptrs to swig ptrs
		set dtPtr [ blobutils_cast_int_to_dbl_ptr $dtPtr]
		set deltaPtr [ blobutils_cast_int_to_dbl_ptr $deltaPtr]
		set gammaPtr [ blobutils_cast_int_to_dbl_ptr $gammaPtr]
		set dimsPtr [ blobutils_cast_int_to_int_ptr $dimsPtr]
		#Load library
		load ../bld/libcohmm_swiftt.so
		puts $doKriging
		#Call function
		set retValue [ prepFirstFlux $doKriging $doCoMD $dimsPtr $dtPtr $deltaPtr $gammaPtr $step ]
		#Return retValue
		return $retValue
	}

	proc prepSecondWrapper {doKriging doCoMD bDims bDt bDelta bGamma step} {
		#Unpack blobs into lists
		set deltaPtr [ lindex $bDelta 0]
		set dimsPtr [ lindex $bDims 0 ]
		set dtPtr [ lindex $bDt 0 ]
		set gammaPtr [ lindex $bGamma 0]
		#Convert ptrs to swig ptrs
		set dtPtr [ blobutils_cast_int_to_dbl_ptr $dtPtr]
		set deltaPtr [ blobutils_cast_int_to_dbl_ptr $deltaPtr]
		set gammaPtr [ blobutils_cast_int_to_dbl_ptr $gammaPtr]
		set dimsPtr [ blobutils_cast_int_to_int_ptr $dimsPtr]
		#Load library
		load ../bld/libcohmm_swiftt.so
		#Call function
		set retValue [ prepSecondFlux $doKriging $doCoMD $dimsPtr $dtPtr $deltaPtr $gammaPtr $step ]
		#Return retValue
		return $retValue
	}

	proc prepThirdWrapper {doKriging doCoMD bDims bDt bDelta bGamma step} {
		#Unpack blobs into lists
		set deltaPtr [ lindex $bDelta 0]
		set dimsPtr [ lindex $bDims 0 ]
		set dtPtr [ lindex $bDt 0 ]
		set gammaPtr [ lindex $bGamma 0]
		#Convert ptrs to swig ptrs
		set dtPtr [ blobutils_cast_int_to_dbl_ptr $dtPtr]
		set deltaPtr [ blobutils_cast_int_to_dbl_ptr $deltaPtr]
		set gammaPtr [ blobutils_cast_int_to_dbl_ptr $gammaPtr]
		set dimsPtr [ blobutils_cast_int_to_int_ptr $dimsPtr]
		#Load library
		load ../bld/libcohmm_swiftt.so
		#Call function
		set retValue [ prepThirdFlux $doKriging $doCoMD $dimsPtr $dtPtr $deltaPtr $gammaPtr $step ]
		#Return retValue
		return $retValue
	}

	proc prepLastWrapper {doKriging doCoMD bDims bDt bDelta bGamma step} {
		#Unpack blobs into lists
		set deltaPtr [ lindex $bDelta 0]
		set dimsPtr [ lindex $bDims 0 ]
		set dtPtr [ lindex $bDt 0 ]
		set gammaPtr [ lindex $bGamma 0]
		#Convert ptrs to swig ptrs
		set dtPtr [ blobutils_cast_int_to_dbl_ptr $dtPtr]
		set deltaPtr [ blobutils_cast_int_to_dbl_ptr $deltaPtr]
		set gammaPtr [ blobutils_cast_int_to_dbl_ptr $gammaPtr]
		set dimsPtr [ blobutils_cast_int_to_int_ptr $dimsPtr]
		#Load library
		load ../bld/libcohmm_swiftt.so
		#Call function
		set retValue [ prepLastFlux $doKriging $doCoMD $dimsPtr $dtPtr $deltaPtr $gammaPtr $step ]
		#Return retValue
		return $retValue
	}

	proc finishStepWrapper {doKriging doCoMD bDims bDt bDelta bGamma step} {
		#Unpack blobs into lists
		set deltaPtr [ lindex $bDelta 0]
		set dimsPtr [ lindex $bDims 0 ]
		set dtPtr [ lindex $bDt 0 ]
		set gammaPtr [ lindex $bGamma 0]
		#Convert ptrs to swig ptrs
		set dtPtr [ blobutils_cast_int_to_dbl_ptr $dtPtr]
		set deltaPtr [ blobutils_cast_int_to_dbl_ptr $deltaPtr]
		set gammaPtr [ blobutils_cast_int_to_dbl_ptr $gammaPtr]
		set dimsPtr [ blobutils_cast_int_to_int_ptr $dimsPtr]
		#Load library
		load ../bld/libcohmm_swiftt.so
		#Call function
		set retValue [ finishStep $doKriging $doCoMD $dimsPtr $dtPtr $deltaPtr $gammaPtr $step ]
		#Return retValue
		return $retValue
	}

	proc prepVTK {doKriging doCoMD bDims bDt bDelta bGamma step} {
		#Unpack blobs into lists
		set deltaPtr [ lindex $bDelta 0]
		set dimsPtr [ lindex $bDims 0 ]
		set dtPtr [ lindex $bDt 0 ]
		set gammaPtr [ lindex $bGamma 0]
		#Convert ptrs to swig ptrs
		set dtPtr [ blobutils_cast_int_to_dbl_ptr $dtPtr]
		set deltaPtr [ blobutils_cast_int_to_dbl_ptr $deltaPtr]
		set gammaPtr [ blobutils_cast_int_to_dbl_ptr $gammaPtr]
		set dimsPtr [ blobutils_cast_int_to_int_ptr $dimsPtr]
		#Load library
		load ../bld/libcohmm_swiftt.so
		#Call function
		set retValue [ outputVTK $doKriging $doCoMD $dimsPtr $dtPtr $deltaPtr $gammaPtr $step ]
		#Return retValue
		return $retValue
	}

	proc cloudFluxWrapper {doKriging doCoMD step phase task} {
		#Load library
		load ../bld/libcohmm_swiftt.so
		#Call function
		set retValue [ cloudFlux $doKriging $doCoMD $step $phase $task ]
		#Return retValue
		return $retValue
	}

	proc shortCircuitWrapper {bDims step} {
		#Unpack blobs into lists
		set dimsPtr [ lindex $bDims 0 ]
		#Convert ptrs to swig ptrs
		set dimsPtr [ blobutils_cast_int_to_int_ptr $dimsPtr]
		#Load library
		load ../bld/libcohmm_swiftt.so
		#Call function
		set retValue [ tryShortCircuit $dimsPtr $step ]
		#Return retValue
		return $retValue
	}

}

package provide cohmm_swiftt $cohmm_swiftt::version

