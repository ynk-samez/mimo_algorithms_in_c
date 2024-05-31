# 実行ファイル
TARGET = main.app

# ソースファイル
INCDIR  = -I../include/
SOURCES = main.c ../include/complexArith.c ../include/def.c ../include/Rfunction.c ../include/error_detection.c
# オブジェクトファイル
OBJECTS = $(SOURCES:.c=.o)

# コンパイラ
CC = gcc
# コンパイルオプション
CFLAGS = -O0 -Wall 

# 実行ファイルのビルド
$(TARGET): $(OBJECTS) $(SOURCES)
	$(CC) $(CFLAGS)  $(INCDIR) -c $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) -o $(TARGET)
	./$(TARGET)

# オブジェクトファイルのビルド
%.o: %.c
	$(CC) $(CFLAGS)  $(INCDIR) -c $< -o $@

# 不要ファイルの削除
clean:
	rm -rf *.gch *.o $(TARGET)


debug:$(SOURCES)
	$(CC) -g $(SOURCES) -o seg
	gdb ./seg

