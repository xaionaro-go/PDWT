package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"net/http"

	"github.com/pierrepaleo/PDWT/src/service/types"
)

func assertNoError(err error) {
	if err != nil {
		panic(err)
	}
}

func main() {
	backendPathFlag := flag.String("backend-path", "./service_backend", "")
	listenAddrFlag := flag.String("listen-addr", ":13407", "")
	flag.Parse()

	backend, err := getBackend(*backendPathFlag)
	assertNoError(err)

	http.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
		body, err := ioutil.ReadAll(r.Body)
		if err != nil {
			w.WriteHeader(400)
			w.Write([]byte(fmt.Sprintf("unable to read body: %v", err)))
			return
		}

		var request types.Request
		err = json.Unmarshal(body, &request)
		if err != nil {
			w.WriteHeader(400)
			w.Write([]byte(fmt.Sprintf("unable to un-JSON-ize: %v", err)))
			return
		}

		response, err := backend.Execute(&request)
		if err != nil {
			w.WriteHeader(500)
			w.Write([]byte(fmt.Sprintf("unable to process the request: %v", err)))
			return
		}

		responseJSON, err := json.Marshal(response)
		if err != nil {
			w.WriteHeader(500)
			w.Write([]byte(fmt.Sprintf("unable to JSON-ize the response: %v", err)))
			return
		}

		w.Header().Set("Content-Type", "application/json")
		w.Write(responseJSON)
	})

	log.Printf("listening for JSON requests at http://%s", *listenAddrFlag)
	err = http.ListenAndServe(*listenAddrFlag, nil)
	assertNoError(err)
}
