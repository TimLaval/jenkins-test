pipeline {
    agent {
        dockerfile true
    }
    stages {
        stage('Git clone') {
                git url: 'https://github.com/TimLaval/jenkins-test'
        }

        stage('Build slave') {
            steps {
                sh 'node --version'    
            }
        }
    }
}
