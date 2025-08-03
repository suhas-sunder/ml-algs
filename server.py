# server.py
import socket

def run_server():
    server_sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_sock.bind(("localhost", 12345))
    server_sock.listen(1)
    print("Server listening on port 12345...")

    client_sock, addr = server_sock.accept()
    print("Connected to", addr)

    while True:
        data = client_sock.recv(1024)
        if not data:
            break
        print("Received:", data.decode())
        client_sock.send(data)

    client_sock.close()
    server_sock.close()

run_server()
