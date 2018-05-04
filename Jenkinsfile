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

		//Retrieve jenkins UID and GID on host node.

		docker.image("${PROPERTIES['DOCKER_IMAGE_NAME']}:latest).inside("-u root --name ${PROPERTIES['DOCKER_IMAGE_NAME']} -v ${PROPERTIES['ABCD_PATH']}:/media -e HOST_USER_ID=${jenkins_uid} -e HOST_USER_GID=${jenkins_gid}"){
			     sh "/runTest.sh"}

}
			
