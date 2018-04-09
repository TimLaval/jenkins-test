node {
    stage "Prepare environment"
          checkout scm
          def env = docker.build 'debian-slave'
    
          env.inside {
              stage("Checkout and build deps"){
                sh "echo 'DEBUG'"
              }
              
              stage("Doxygen"){
                publishHTML target: [$class: 'HtmlPublisherTarget', reportName: 'Doxygen', reportDir: 'build/linux/doc/html', reportFiles: 'index.html']
              }
           
          }
}
