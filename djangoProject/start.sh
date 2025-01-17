#!/bin/bash

screen -dmS my_session sshpass -p 'Zyf128704' ssh -o "ServerAliveInterval 5" -o "ServerAliveCountMax 10" -R 8324:10.40.16.247:8000 -g root@47.92.111.142

# 激活虚拟环境
conda activate DjangoprojectUnregistered

# 运行后端保持在后台
screen -dmS my_session1 daphne -b 10.40.16.247 -p 8000 djangoProject.asgi:application


