node {
    stage "Prepare environment"
          checkout scm
          def env = docker.build 'debian slave'
    
          env.inside {
            stage "Checkout and build deps"
                sh "echo 'DEBUG'"

           
          }
}
