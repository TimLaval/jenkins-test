### Pipeline jenkins to execute container 

xhost +local:root
sudo docker run -it --rm --net=host -e DISPLAY=$DISPLAY -v $(pwd)/:/work/abcd/dermscan -v /tmp/.X11-unix:/tmp/.X11-unix abcd:latest

