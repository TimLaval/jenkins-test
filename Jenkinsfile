node('master') { 
        stage('Git clone') {
            git url: 'https://github.com/TimLaval/jenkins-test'
        }
        
        stage('Debian build') {
            app = docker.build("agent")
        }
        
        stage("Test image"){
            app.inside {
                sh 'echo "Test passed"'
            }
        }
    
}
