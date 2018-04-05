node('master'){
    def app 
    agent {
        dockerfile true
    }
  
    stage('Git clone') {
            git url: 'https://github.com/TimLaval/jenkins-test'
    }

    stage('Build slave') {
        steps {
            sh 'node --version'
        }
    }
}
