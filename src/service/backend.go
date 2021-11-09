package main

import "C"
import (
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"os/exec"
	"reflect"
	"strconv"
	"sync"
	"unsafe"

	"github.com/pierrepaleo/PDWT/src/service/types"
)

type backend struct {
	sync.Mutex
	*exec.Cmd
	BackendExecutablePath string
	Stdin                 io.WriteCloser
	Stdout                io.ReadCloser
	Stderr                io.ReadCloser
}

func getBackend(backendExecutablePath string) (*backend, error) {
	b := &backend{
		BackendExecutablePath: backendExecutablePath,
	}
	if err := b.start(); err != nil {
		return nil, fmt.Errorf("unable to start: %w", err)
	}
	return b, nil
}

func (b *backend) start() error {
	cmd := exec.Command(b.BackendExecutablePath)
	stdin, err := cmd.StdinPipe()
	if err != nil {
		return fmt.Errorf("unable to get stdin pipe: %w", err)
	}
	stdout, err := cmd.StdoutPipe()
	if err != nil {
		return fmt.Errorf("unable to get stdin pipe: %w", err)
	}
	stderr, err := cmd.StderrPipe()
	if err != nil {
		return fmt.Errorf("unable to get stderr pipe: %w", err)
	}

	err = cmd.Start()
	if err != nil {
		return fmt.Errorf("unable to start a backend: %w", err)
	}

	b.Cmd = cmd
	b.Stdin = stdin
	b.Stdout = stdout
	b.Stderr = stderr
	return nil
}

func (b *backend) Close() {
	b.Mutex.Lock()
	defer b.Mutex.Unlock()

	b.close()
}

func (b *backend) close() {
	b.Stdin.Close()
	b.Stdout.Close()
	b.Process.Kill()
}

func (b *backend) restart() {
	b.close()
	b.start()
}

func (b *backend) Execute(req *types.Request) (*types.Response, error) {
	tmpDir := os.TempDir()

	requestFile, err := ioutil.TempFile(tmpDir, ".PDWT-request-")
	if err != nil {
		return nil, fmt.Errorf("unable to create a request file: %w")
	}

	responseFile, err := ioutil.TempFile(tmpDir, ".PDWT-response-")
	if err != nil {
		return nil, fmt.Errorf("unable to create a response file: %w")
	}

	b.Lock()
	defer b.Unlock()

	valuesSlice := (*reflect.SliceHeader)(unsafe.Pointer(&req.Values))
	valuesSlice.Cap *= int(unsafe.Sizeof(float32(0)))
	valuesSlice.Len *= int(unsafe.Sizeof(float32(0)))
	n, err := requestFile.Write(*(*[]byte)((unsafe.Pointer)(valuesSlice)))
	if err != nil {
		return nil, fmt.Errorf("unable to write data to the request file: %w", err)
	}
	if n != valuesSlice.Len {
		return nil, fmt.Errorf("wrote invalid length of data to the request file: %d != %d", n, valuesSlice.Len)
	}

	const (
		fieldSeparator    = '\000'
		requestTerminator = '\n'
	)

	var buf bytes.Buffer
	_, err = buf.WriteString(req.Wavelet)
	assertNoError(err)
	err = buf.WriteByte(fieldSeparator)
	assertNoError(err)
	_, err = buf.WriteString(strconv.FormatUint(uint64(req.Levels), 10))
	assertNoError(err)
	err = buf.WriteByte(fieldSeparator)
	assertNoError(err)
	_, err = buf.WriteString(strconv.FormatBool(req.IsInverse))
	assertNoError(err)
	err = buf.WriteByte(fieldSeparator)
	assertNoError(err)
	_, err = buf.WriteString(requestFile.Name())
	assertNoError(err)
	err = buf.WriteByte(fieldSeparator)
	assertNoError(err)
	_, err = buf.WriteString(responseFile.Name())
	err = buf.WriteByte(requestTerminator)
	assertNoError(err)

	_, err = b.Stdin.Write(buf.Bytes())
	if err != nil {
		b.restart()
		return nil, fmt.Errorf("unable to write request '%s' to the backend's stdin: %w", buf.String(), err)
	}

	var outcome [256]byte
	n, err = b.Stdout.Read(outcome[:])
	if err != nil {
		var errorDescriptionBuf [1024]byte
		var errorDescription []byte
		n, _ := b.Stderr.Read(errorDescriptionBuf[:])
		if n > 0 {
			errorDescription = errorDescriptionBuf[:n]
		}
		b.restart()
		return nil, fmt.Errorf("unable to read the outcome from the backend's stdout: %w, stderr: %s", err, errorDescription)
	}

	if outcome[0] != '\n' {
		b.restart()
		return nil, fmt.Errorf("unsupported backend protocol, received: %X:%s", outcome[:n], outcome[:n])
	}

	responseRawData, err := ioutil.ReadAll(responseFile)
	if err != nil {
		b.restart()
		return nil, fmt.Errorf("read raw data from the response file: %w", err)
	}

	if len(responseRawData)%int(unsafe.Sizeof(float32(0))) != 0 {
		b.restart()
		return nil, fmt.Errorf("unexpected size: %d %% 4 != 0", len(responseRawData))
	}

	responseSlice := (*reflect.SliceHeader)(unsafe.Pointer(&responseRawData))
	responseSlice.Len /= int(unsafe.Sizeof(float32(0)))
	responseSlice.Cap /= int(unsafe.Sizeof(float32(0)))
	responseData := *(*[]float32)(unsafe.Pointer(responseSlice))

	return &types.Response{
		Values: responseData,
	}, nil
}
