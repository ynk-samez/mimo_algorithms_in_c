# 実行ファイル
TARGET = main.app

# ソースファイル
SOURCES = main.c complexArith.c def.c Rfunction.c
# オブジェクトファイル
OBJECTS = $(SOURCES:.c=.o)

# コンパイラ
CC = gcc
# コンパイルオプション
CFLAGS = -O0 -Wall 

# 実行ファイルのビルド
$(TARGET): $(OBJECTS) $(SOURCES)
	$(CC) $(CFLAGS) -c $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) -o $(TARGET)
	./$(TARGET)

# オブジェクトファイルのビルド
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# 不要ファイルの削除
clean:
	rm -rf *.gch *.o $(TARGET)


debug:$(SOURCES)
	$(CC) -g $(SOURCES) -o seg
	gdb ./seg
