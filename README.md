# SimpleRayTracing

本项目基于 Kevin Beason 的 [smallpt开源项目](http://www.kevinbeason.com/smallpt/)，以纯CPU渲染的方式实现路径追踪。

# 编译方法

`g++ main.cpp -O3 -fopenmp -o main`

# 渲染结果演示：

Sample1: 1024 * 768 @ 5000samps

`./main 5000`

![Sample1.jpg](./Sample1.jpg)

Sample2: 1920 * 1080 @ 25000samps

`./main 25000`

![Sample2.jpg](./Sample2.jpg)



