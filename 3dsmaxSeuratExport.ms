----------------------------------------------------------------------------------------------------------------
--
--		3ds Max Seurat Export 
--		Generates a JSON manifest and 3ds Max camera rig for Seurat.
--
--		Based on "seurat_rig.py" from the GoogleVR Seurat Github distro.
--
--		(c) SuperRune 2018
--		rune@superrune.com
--
----------------------------------------------------------------------------------------------------------------

struct seuratExport
(
	prefix = "seurat_cam",
	faceNames = #("front", "back", "left", "right", "bottom", "top"),
	depthChannel = "Z.R",
	depthType = "EYE_Z",
	
	function verbose verboseText =
	(
		format "-- %\n" verboseText
		return True
	),

	-- Orientation matrices for the face cameras
	function worldFromEyeMatrixFromFace faceName inPos =
	(
		if faceName == "front" then (
			return #( 1.0,  0.0,  0.0,  inPos[1], \
					  0.0,  1.0,  0.0,  inPos[2], \
					  0.0,  0.0,  1.0,  inPos[3], \
					  0.0,  0.0,  0.0,  1.0 )
		) else if faceName == "back" then (
			return #(-1.0,  0.0,  0.0,  inPos[1], \
					  0.0,  1.0,  0.0,  inPos[2], \
					  0.0,  0.0, -1.0,  inPos[3], \
					  0.0,  0.0,  0.0,  1.0 )
		) else if faceName == "left" then (
			return #( 0.0,  0.0,  1.0,  inPos[1], \
					  0.0,  1.0,  0.0,  inPos[2], \
					 -1.0,  0.0,  0.0,  inPos[3], \
					  0.0,  0.0,  0.0,  1.0 )
		) else if faceName == "right" then (
			return #( 0.0,  0.0, -1.0,  inPos[1], \
					  0.0,  1.0,  0.0,  inPos[2], \
					  1.0,  0.0,  0.0,  inPos[3], \
					  0.0,  0.0,  0.0,  1.0 )
		) else if faceName == "bottom" then (
			return #( 1.0,  0.0,  0.0,  inPos[1], \
					  0.0,  0.0,  1.0,  inPos[2], \
					  0.0, -1.0,  0.0,  inPos[3], \
					  0.0,  0.0,  0.0,  1.0 )
		) else if faceName == "top" then (
			return #( 1.0,  0.0,  0.0,  inPos[1], \
					  0.0,  0.0, -1.0,  inPos[2], \
					  0.0,  1.0,  0.0,  inPos[3], \
					  0.0,  0.0,  0.0,  1.0 )
		) else return False
	),

	function pointInRadius inPos radius sample =
	(
		zeroPos = [inPos.x-radius, inPos.y-radius, inPos.z-radius]
		newPos = zeroPos + [radius*2*sample[1], radius*2*sample[2], radius*2*sample[3]]
		return newPos
	),
	
	-- Creates a cube-face 90 degree FOV OpenGL-style projection matrix
	function cubeFaceProjectionMatrix_old near far =
	(
		near = near as float
		far = far as float
		
		left = -near
		right = near
		bottom = -near
		top = near
		a = (2.0 * near) / (right - left)
		b = (2.0 * near) / (top - bottom)
		c = (right + left) / (right - left)
		d = (top + bottom) / (top - bottom)
		ee = (near + far) / (near - far)
		f = (2.0 * near * far) / (near - far)

		return #(	a,   0.0,     c,  0.0, \
					0.0,    b,    d,  0.0, \
					0.0,  0.0,   ee,    f, \
					0.0,  0.0, -1.0,  0.0 )
	),
	
	-- Creates a cube-face 90 degree FOV OpenGL-style projection matrix, shorter version
	function cubeFaceProjectionMatrix near far =
	(
		far = far as float
		near = near as float
		
		a = (-(far + near)) / (far-near)
		b = (-2.0*far*near) / (far-near)
		return #(	1.0,  0.0,  0.0,  0.0, \
					0.0,  1.0,  0.0,  0.0, \
					0.0,  0.0,    a,    b, \
					0.0,  0.0, -1.0,  0.0 )
	),
	
	-- Calculates the radical inverse of |a| in base |base|
	function radicalInverse a base =
	(
		reversed_digits = 0
		base_n = 1
		
		while a > 0 do
		(
			next_a = a / base
			digit = a - next_a * base
			reversed_digits = reversed_digits * base + digit
			base_n = base_n * base
			a = next_a
		)
		
		return amin #((reversed_digits / base_n as float), 1.0)
	),

	
	-- Generates camera positions scattered in a radius around a position
	function generateCameraPositions inPos numPos radius = 
	(
		samples = #()
		max_sample = [0.0, 0.0, 0.0]
		
		for c = 0 to numPos-1 do
		(
			-- Use a 3D Hammersley point set for the samples
			sample = #(c / numPos as float, radicalInverse c 2, radicalInverse c 3)
			for dim = 1 to 3 do
			(
				max_sample[dim] = amax #(max_sample[dim], sample[dim])
			)
			append samples sample
		)
		
		closestSample = numPos
		closestDist = 2000
		
		for c = 1 to numPos do
		(
			curDist = distance [0.5,0.5,0.5] [samples[c][1], samples[c][2], samples[c][3]]
			if curDist < closestDist then
			(
				closestDist = curDist
				closestSample = c
			)
		)
		
		-- Remove the sample with the position closest to center, add a center position in the first array entry
		verbose ("Camera sample " + closestSample as string + ":" + samples[closestSample] as string + " is closest to headBox center.")
		deleteItem samples closestSample 
		insertItem [0.5,0.5,0.5] samples 1
		
		camera_positions = #()
		for sample in samples do append camera_positions (pointInRadius inPos radius sample)
		
		return camera_positions
	),
	
	function createCameras cameraPositions nearClip:1.0 = 
	(
		for c = 1 to cameraPositions.count do
		(
			camPos = (cameraPositions[c])
			camParent = point name:(prefix + (formattedPrint c format:"03d")) pos:camPos size:1.0
			for face in faceNames do 
			(
				faceCam = Freecamera fov:90 nearClip:nearClip farClip:farClip name:(camParent.name + "_" + face) pos:camPos showCone:true
				faceCam.parent = camParent
				if face == "front" then rotate faceCam (angleaxis 90 [1,0,0])
				else if face == "back" then (rotate faceCam (angleaxis 90 [1,0,0]); rotate faceCam (angleaxis 180 [0,0,1]))
				else if face == "left" then (rotate faceCam (angleaxis 90 [1,0,0]); rotate faceCam (angleaxis 90 [0,0,1]))
				else if face == "right" then (rotate faceCam (angleaxis 90 [1,0,0]); rotate faceCam (angleaxis -90 [0,0,1]))
				else if face == "top" then rotate faceCam (angleaxis 180 [1,0,0])
			)
		)
		return True
	),
	
	function createAnimatedCamera cameraPositions frameOffset:0 nearClip:1.0 farClip:1000000.0 =
	(
		faceCam = Freecamera fov:90 nearClip:nearClip farClip:farClip name:(prefix + "_render") pos:[0,0,0] showCone:true
		
		for c = 1 to cameraPositions.count do point name:(prefix + "_" + (formattedPrint (c-1) format:"03d")) pos:(cameraPositions[c]) size:1.0 centermarker:true axistripod:false cross:false box:false
		
		-- Set frame range
		animationRange = interval (frameOffset) ((cameraPositions.count*6)+frameOffset-1)
		
		with animate on
		(
			-- Animate rotation of camera to cover all faces of the headbox
			sliderTime = 0+frameOffset
			rotate faceCam (angleaxis 90 [1,0,0]) --Front
			sliderTime = 1+frameOffset
			rotate faceCam (angleaxis 180 [0,0,1]) --Back
			sliderTime = 2+frameOffset
			rotate faceCam (angleaxis -90 [0,0,1]) --Left
			sliderTime = 3+frameOffset
			rotate faceCam (angleaxis 180 [0,0,1]) --Right
			sliderTime = 4+frameOffset
			rotate faceCam (angleaxis 90 [0,0,1])
			rotate faceCam (angleaxis -90 [1,0,0]) --Top
			sliderTime = 5+frameOffset
			rotate faceCam (angleaxis 180 [1,0,0]) --Bottom
			sliderTime = 6+frameOffset
			rotate faceCam (angleaxis -90 [1,0,0])
			--rotate faceCam (angleaxis 180 [0,0,1])
			
			-- Animate position of camera to cover all generated positions
			for t = 1 to cameraPositions.count do
			(
				sliderTime = t+((t-1)*5)+frameOffset-1
				faceCam.pos = cameraPositions[t]
			)
		)
		
		-- Non-default timelines create an extra key that messes up looping, so we delete it
		if frameOffset != 0 then deleteKey faceCam.rotation.controller.keys 1
		
		-- Make rotations repeat and positions step
		setAfterORT faceCam.rotation.controller.X_Rotation.controller #cycle
		setAfterORT faceCam.rotation.controller.Y_Rotation.controller #cycle
		setAfterORT faceCam.rotation.controller.Z_Rotation.controller #cycle
		for k in faceCam.position.controller.X_Position.controller.keys do k.inTangentType = k.outTangentType = #step
		for k in faceCam.position.controller.Y_Position.controller.keys do k.inTangentType = k.outTangentType = #step
		for k in faceCam.position.controller.Z_Position.controller.keys do k.inTangentType = k.outTangentType = #step

		sliderTime = frameOffset
		select faceCam
		return faceCam
	),
	
	function formatMatrixForJSON matrixArray numTabs:7 =
	(
		matrixString = ""
		tabString = ""
		for c = 1 to numTabs do tabString += "\t"
		for c = 1 to 4 do
		(
			offset = (c-1) * 4
			matrixString += tabString + matrixArray[offset+1] as string + ", " + \
				matrixArray[offset+2] as string + ", " + matrixArray[offset+3] as string + \
				", " + matrixArray[offset+4] as string 
			if c != 4 then matrixString += ",\n" else matrixString += "\n"
		)
		
		return matrixString
	),
	
	function writeJSON cameraPositions headBoxCenter targetFolder frameOffset:0 nearClip:1.0 farClip:1000000.0 resolution:1024 numPos:32 =
	(
		outFile = createFile (targetFolder + "manifest.json" )
		frameCounter = frameOffset
		
		-- No JSON support built into maxscript. Now cover your eyes, this is going to get dirty.
		format "{\n" to:outFile
		format "\t\"settingsinfo\": [ \"Exported from 3ds Max. Headbox center is %, near clip is %, far clip is: %.\"],\n" (headBoxCenter as string) (nearClip as string) (farClip as string) to:outFile
		format "\t\"view_groups\": [\n\t\t{\n" to:outFile
		for c = 1 to cameraPositions.count do
		(
			-- Positions are relative to headbox center
			position = [cameraPositions[c][1] - headBoxCenter[1], cameraPositions[c][2] - headBoxCenter[2], cameraPositions[c][3] - headBoxCenter[3]]
			-- Convert to left-handed coordinate system
			position = [position[1], position[3], -position[2]]
				
			format "\t\t\t\"views\": [\n" to:outFile
			
			for face in faceNames do
			(
				-- Camera matrices
				format "\t\t\t\t{\n\t\t\t\t\t\"projective_camera\": {\n" to:outFile
				format "\t\t\t\t\t\t\"image_width\": %,\n" (resolution as string) to:outFile
				format "\t\t\t\t\t\t\"image_height\": %,\n" (resolution as string) to:outFile
				format "\t\t\t\t\t\t\"clip_from_eye_matrix\": [\n" to:outFile
				format "%" (formatMatrixForJSON (cubeFaceProjectionMatrix nearClip farClip) numTabs:7) to:outFile
				format "\t\t\t\t\t\t],\n" to:outFile
				format "\t\t\t\t\t\t\"world_from_eye_matrix\": [\n" to:outFile
				format "%" (formatMatrixForJSON (worldFromEyeMatrixFromFace face position) numTabs:7) to:outFile
				format "\t\t\t\t\t\t],\n" to:outFile
				format "\t\t\t\t\t\t\"depth_type\": \"%\"\n" depthType to:outFile
				format "\t\t\t\t\t},\n" to:outFile
				
				-- File locations
				format "\t\t\t\t\t\"depth_image_file\": {\n" to:outFile
				format "\t\t\t\t\t\t\"color\": {\n" to:outFile
				format "\t\t\t\t\t\t\t\"path\": \"%\",\n" (prefix + "." + (formattedPrint frameCounter format:"04d") + ".exr") to:outFile
				format "\t\t\t\t\t\t\t\"channel_0\": \"R\",\n" to:outFile
				format "\t\t\t\t\t\t\t\"channel_1\": \"G\",\n" to:outFile
				format "\t\t\t\t\t\t\t\"channel_2\": \"B\",\n" to:outFile
				format "\t\t\t\t\t\t\t\"channel_alpha\": \"CONSTANT_ONE\"\n" to:outFile
				format "\t\t\t\t\t\t},\n" to:outFile
				format "\t\t\t\t\t\t\"depth\": {\n" to:outFile
				format "\t\t\t\t\t\t\t\"path\": \"%\",\n" (prefix + "." + (formattedPrint frameCounter format:"04d") + ".exr") to:outFile
				format "\t\t\t\t\t\t\t\"channel_0\": \"%\"\n" depthChannel to:outFile
				format "\t\t\t\t\t\t}\n\t\t\t\t\t}\n" to:outFile
				
				if face == faceNames[6] then
					format "\t\t\t\t}\n" to:outFile
				else
					format "\t\t\t\t},\n" to:outFile
				
				frameCounter += 1
			
			)
			if c != cameraPositions.count then format "\t\t\t]\n\t\t},\n\t\t{\n" to:outFile
			
		)
		format "\t\t\t]\n\t\t}\n\t]\n}" to:outFile
		
		close outFile
		
		return True
	),
	
	function main headBoxCenter headBoxSize targetFolder frameOffset:0 nearClip:1.0 farClip:10000.0 resolution:1024 numPos:32 =
	(
		-- Set render resolution and render out name
		renderWidth = renderHeight = resolution
		rendOutputFilename = targetFolder + prefix + ".exr"
		
		-- Generate all cameras and output JSON
		cameraPositions = generateCameraPositions headBoxCenter numPos headBoxSize
		createAnimatedCamera cameraPositions frameOffset:frameOffset
		writeJSON cameraPositions headBoxCenter targetFolder frameOffset:frameOffset nearClip:nearClip farClip:farClip resolution:resolution numPos:numPos 
		
		return True
	)
)

seuratExport = seuratExport()

rollout seuratExport_UI "Seurat Exporter"
(
	group ""
	(
		label lbl0 "Seurat Exporter - Superrune 2018" height:24 offset:[0,12]
	)
	
	editText nameTxt "Name: " fieldWidth:150 text:(substituteString (getFilenameFile maxFileName) " " "") offset:[0,7] align:#center
	editText pathTxt "Path: " fieldWidth:150 text:(maxFilePath + "Output\\") offset:[58,3] align:#center across:2
	button pathLoad "..." width:16 offset:[50,1]
	
	group ""
	(
		button getPos "Get position from selection" width:180
		spinner camX 	"X:" 	range:[-50000,50000,0] fieldWidth:40 offset:[0,0] align:#center across:3
		spinner camY 	"Y:" 	range:[-50000,50000,0] fieldWidth:40 offset:[0,0] align:#center
		spinner camZ 	"Z:" 	range:[-50000,50000,50] fieldWidth:40 offset:[0,0] align:#center
		
		spinner viewBoxRad "Viewbox Radius       "  range:[0,90000,10] fieldWidth:50 align:#center offset:[0,10]
	)
	
	group ""
	(
		label lbl2 "Render Resolution" width:100 across:2 offset:[30,0]
		dropDownList resolution "" items:#("512","1024","2048","8192") selection:2 width:60 offset:[18, -3] 
		label lbl3 "Camera Positions" width:100 across:2 offset:[30,4]
		dropDownList camPos "" items:#("8","16","32","64") selection:3 width:60 offset:[18, 1]  
	)
	
	group ""
	(
		spinner frameOffset "Frame Start            " range:[0,90000,animationRange.start] fieldWidth:50 type:#integer align:#center
		spinner clipStart "Camera clip" range:[0,100000000000,1] fieldWidth:50 type:#float align:#center offset:[30,0] across:2
		spinner clipEnd "" range:[0,100000000000,10000] fieldWidth:50 type:#float align:#center offset:[10,0] 
	)
	
	button doExport "Generate" width:180 offset:[0,10] height:40
	
	-- UI handlers
	on seuratExport_UI open do 
	(
		if maxFileName == "" then nameTxt.text = "Viewbox"
		seuratExport.prefix = nameTxt.text + "_seurat_cam"
		if maxFilePath == "" then pathTxt.text = "C:\\Temp\\Output\\"
	)
	
	on pathLoad pressed do
	(
		local newDir = getSavePath caption:"Select an output path." initialDir:pathTxt.text
		if newDir != undefined then pathTxt.text = newDir + "\\"
	)
		
	on nameTxt changed txt do 
	(
		nameTxt.text = substituteString txt " " ""
		seuratExport.prefix = txt + "_seurat_cam"
	)
	
	on getPos pressed do 
	(
		print ("Setting " + [camX.value, camY.value, camZ.value] as string + " to " + $.center as string)
		camX.value = $.center[1]
		camY.value = $.center[2]
		camZ.value = $.center[3]
	)
	
	on doExport pressed do 
	(
		seuratExport.verbose ("Creating " + camPos.selected + " camera positions and JSON for '" + nameTxt.text + "' to path '" + pathTxt.text + "'.")
		if (doesDirectoryExist pathTxt.text) == false then
		(
			shouldWeCreate = queryBox ("Path '" + pathTxt.text + "' does not exist. Create it?") title:"Path not found!"
			if shouldWeCreate == True then makeDir pathTxt.text
			else return false
		)
		delete $'*_seurat_cam_*'
		seuratExport.main [camX.value, camY.value, camZ.value] viewBoxRad.value pathTxt.text nearClip:clipStart.value farClip:clipEnd.value \
		frameOffset:frameOffset.value resolution:(resolution.selected as integer) numPos:(camPos.selected as integer)
	)
)

createDialog seuratExport_UI 250 430 pos:[100,100] lockHeight:true