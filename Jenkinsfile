#!groovy

node {

	try {
		checkout()

		buildAndTest()

	}catch(Exception err) {
		echo "Something went wrong: ${err}"
	}
}

def checkout() {
	stage('Checkout code'){
		checkout scm
		PROPERTIES = readProperties file: 'jenkins.properties'
	}
}

def buildAndTest() {
	stage('Building and Testing') {
		// Retrieve current workspace directory.
		PWD = pwd()

		// Retrieve jenkins UID and GID on host node.
		jenkins_uid = sh(returnStdout: true, script: "cat /etc/passwd | grep 'jenkins' | cut -d : -f 3").trim()
		jenkins_gid = sh(returnStdout: true, script: "cat /etc/passwd | grep 'jenkins' | cut -d : -f 4").trim()


		// Run docker image (+ perform host independent UID/GID mapping).
		docker.image("${PROPERTIES['DOCKER_IMAGE_NAME']}:latest").inside("-u root --name ${PROPERTIES['DOCKER_IMAGE_NAME']} -v ${PROPERTIES['LUMOS_AOSP_PATH']}:/media -e HOST_USER_ID=${jenkins_uid} -e HOST_USER_GID=${jenkins_gid}") {

			// Host independent UID/GID mapping.
			sh "/runTest.sh"

		}
	}

}
			
