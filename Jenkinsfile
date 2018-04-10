node {
    stage "Prepare environment"
          checkout scm
          def env = docker.build 'debian-slave'
    
          env.inside {
              stage("Checkout and build deps"){
                  sh "echo 'DEBUG'"
                  sh 'cmake . -G"Unix Makefiles"'
              }
              
              stage("Doxygen"){
                  sh "make doc"
                  sh "publishHTML target: [$class: 'HtmlPublisherTarget', reportName: 'Doxygen', reportDir: 'build/linux/doc/html', reportFiles: 'index.html']"
              }
           
          }
}
