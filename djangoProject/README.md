1、 执行 pip install -r requirements.txt
可能会存在疏漏的包，在执行时会提示报错没有相应的包，
【解决办法】复制提示的包名，pip install 报名

2、执行sql文件。

3、【启动】
在虚拟环境中的工程目录下 执行 daphne -b 192.168.111.129 -p 8000 djangoProject.asgi:application
【参数解释】要根据你实际的ip的变化
daphne -b 192.168.111.129【IP地址】 -p 8000【端口】 djangoProject.asgi:application