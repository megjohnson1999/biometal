//! Integration tests for HTTP streaming with mock server
//!
//! These tests validate real HTTP behavior using wiremock mock server.

#[cfg(feature = "network")]
mod tests {
    use biometal::io::network::{HttpClient, HttpReader};
    use biometal::BiometalError;
    use std::io::Read;
    use wiremock::matchers::{header, method, path};
    use wiremock::{Mock, MockServer, ResponseTemplate};

    #[tokio::test]
    async fn test_successful_range_request() {
        let mock_server = MockServer::start().await;
        let test_data = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

        // Mock HEAD request
        Mock::given(method("HEAD"))
            .and(path("/test.txt"))
            .respond_with(
                ResponseTemplate::new(200)
                    .insert_header("Content-Length", test_data.len().to_string()),
            )
            .mount(&mock_server)
            .await;

        // Mock range request
        Mock::given(method("GET"))
            .and(path("/test.txt"))
            .and(header("Range", "bytes=0-9"))
            .respond_with(
                ResponseTemplate::new(206)
                    .set_body_bytes(&test_data[0..10])
                    .insert_header("Content-Range", "bytes 0-9/36"),
            )
            .mount(&mock_server)
            .await;

        let url = format!("{}/test.txt", mock_server.uri());

        // Use spawn_blocking to avoid runtime conflicts
        let result = tokio::task::spawn_blocking(move || {
            let client = HttpClient::new()?;
            client.fetch_range(&url, 0, 10)
        })
        .await
        .unwrap();

        let data = result.unwrap();
        assert_eq!(data.as_ref(), &test_data[0..10]);
    }

    #[tokio::test]
    async fn test_server_without_range_support() {
        let mock_server = MockServer::start().await;
        let test_data = b"Full file content";

        // Mock HEAD request
        Mock::given(method("HEAD"))
            .and(path("/no-range.txt"))
            .respond_with(
                ResponseTemplate::new(200)
                    .insert_header("Content-Length", test_data.len().to_string()),
            )
            .mount(&mock_server)
            .await;

        // Server ignores Range header and returns 200
        Mock::given(method("GET"))
            .and(path("/no-range.txt"))
            .respond_with(ResponseTemplate::new(200).set_body_bytes(test_data))
            .mount(&mock_server)
            .await;

        let url = format!("{}/no-range.txt", mock_server.uri());

        let result = tokio::task::spawn_blocking(move || {
            let client = HttpClient::new()?;
            client.fetch_range(&url, 0, 10)
        })
        .await
        .unwrap();

        // Should fail with error about range support
        assert!(result.is_err());
        match result.unwrap_err() {
            BiometalError::Network(msg) => {
                assert!(msg.contains("does not support range requests"));
            }
            e => panic!("Expected Network error, got: {:?}", e),
        }
    }

    #[tokio::test]
    async fn test_range_not_satisfiable() {
        let mock_server = MockServer::start().await;
        let test_data = b"Short";

        // Mock HEAD request
        Mock::given(method("HEAD"))
            .and(path("/short.txt"))
            .respond_with(
                ResponseTemplate::new(200)
                    .insert_header("Content-Length", test_data.len().to_string()),
            )
            .mount(&mock_server)
            .await;

        // Request range beyond file size
        Mock::given(method("GET"))
            .and(path("/short.txt"))
            .respond_with(ResponseTemplate::new(416))
            .mount(&mock_server)
            .await;

        let url = format!("{}/short.txt", mock_server.uri());

        let result = tokio::task::spawn_blocking(move || {
            let client = HttpClient::new()?;
            client.fetch_range(&url, 100, 200)
        })
        .await
        .unwrap();

        // Should fail with range error
        assert!(result.is_err());
        match result.unwrap_err() {
            BiometalError::Network(msg) => {
                assert!(msg.contains("out of bounds"));
            }
            e => panic!("Expected Network error, got: {:?}", e),
        }
    }

    #[tokio::test]
    async fn test_cache_hit() {
        let mock_server = MockServer::start().await;
        let test_data = b"Cached data";

        // Mock HEAD request
        Mock::given(method("HEAD"))
            .and(path("/cached.txt"))
            .respond_with(
                ResponseTemplate::new(200)
                    .insert_header("Content-Length", test_data.len().to_string()),
            )
            .mount(&mock_server)
            .await;

        // Mock range request (should only be called once)
        Mock::given(method("GET"))
            .and(path("/cached.txt"))
            .and(header("Range", "bytes=0-10"))
            .respond_with(
                ResponseTemplate::new(206)
                    .set_body_bytes(&test_data[0..11])
                    .insert_header("Content-Range", "bytes 0-10/11"),
            )
            .expect(1) // Should only be called once
            .mount(&mock_server)
            .await;

        let url = format!("{}/cached.txt", mock_server.uri());

        let (data1, data2, stats) = tokio::task::spawn_blocking(move || {
            let client = HttpClient::new()?;

            // First call - cache miss
            let data1 = client.fetch_range(&url, 0, 11)?;

            // Second call - cache hit
            let data2 = client.fetch_range(&url, 0, 11)?;

            // Get stats
            let stats = client.cache_stats()?;

            Ok::<_, BiometalError>((data1, data2, stats))
        })
        .await
        .unwrap()
        .unwrap();

        assert_eq!(data1.as_ref(), &test_data[0..11]);
        assert_eq!(data2.as_ref(), &test_data[0..11]);
        assert_eq!(stats.entries, 1);
        assert_eq!(stats.current_bytes, 11);
    }

    #[tokio::test]
    async fn test_http_reader_streaming() {
        let mock_server = MockServer::start().await;
        let test_data = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

        // Mock HEAD request
        Mock::given(method("HEAD"))
            .and(path("/stream.txt"))
            .respond_with(
                ResponseTemplate::new(200)
                    .insert_header("Content-Length", test_data.len().to_string()),
            )
            .mount(&mock_server)
            .await;

        // Mock range requests
        Mock::given(method("GET"))
            .and(path("/stream.txt"))
            .respond_with(move |req: &wiremock::Request| {
                let range = req.headers.get("Range").unwrap().to_str().unwrap();
                let range_str = range.strip_prefix("bytes=").unwrap();
                let parts: Vec<&str> = range_str.split('-').collect();
                let start: usize = parts[0].parse().unwrap();
                let end: usize = parts[1].parse().unwrap();

                let chunk = &test_data[start..=end];
                ResponseTemplate::new(206)
                    .set_body_bytes(chunk)
                    .insert_header(
                        "Content-Range",
                        format!("bytes {}-{}/{}", start, end, test_data.len()),
                    )
            })
            .mount(&mock_server)
            .await;

        let url = format!("{}/stream.txt", mock_server.uri());

        let result = tokio::task::spawn_blocking(move || {
            let mut reader = HttpReader::new(&url)?;
            let mut buffer = vec![0u8; 10];
            let mut result = Vec::new();

            loop {
                let n = reader.read(&mut buffer)?;
                if n == 0 {
                    break;
                }
                result.extend_from_slice(&buffer[..n]);
            }

            Ok::<_, BiometalError>(result)
        })
        .await
        .unwrap()
        .unwrap();

        assert_eq!(result.as_slice(), test_data);
    }

    #[tokio::test]
    async fn test_http_reader_eof_detection() {
        let mock_server = MockServer::start().await;
        let test_data = b"Hello";

        // Mock HEAD request
        Mock::given(method("HEAD"))
            .and(path("/eof.txt"))
            .respond_with(
                ResponseTemplate::new(200)
                    .insert_header("Content-Length", test_data.len().to_string()),
            )
            .mount(&mock_server)
            .await;

        // Mock range request
        Mock::given(method("GET"))
            .and(path("/eof.txt"))
            .respond_with(move |req: &wiremock::Request| {
                let range = req.headers.get("Range").unwrap().to_str().unwrap();
                let range_str = range.strip_prefix("bytes=").unwrap();
                let parts: Vec<&str> = range_str.split('-').collect();
                let start: usize = parts[0].parse().unwrap();
                let end: usize = parts[1].parse().unwrap();

                if start >= test_data.len() {
                    ResponseTemplate::new(416)
                } else {
                    let actual_end = end.min(test_data.len() - 1);
                    let chunk = &test_data[start..=actual_end];
                    ResponseTemplate::new(206)
                        .set_body_bytes(chunk)
                        .insert_header(
                            "Content-Range",
                            format!("bytes {}-{}/{}", start, actual_end, test_data.len()),
                        )
                }
            })
            .mount(&mock_server)
            .await;

        let url = format!("{}/eof.txt", mock_server.uri());

        let (n1, n2, n3) = tokio::task::spawn_blocking(move || {
            let mut reader = HttpReader::new(&url)?;
            let mut buffer = vec![0u8; 1024];

            // Read entire file
            let n1 = reader.read(&mut buffer)?;
            assert_eq!(&buffer[..n1], test_data);

            // Next read should return 0 (EOF)
            let n2 = reader.read(&mut buffer)?;

            // Multiple reads at EOF should continue returning 0
            let n3 = reader.read(&mut buffer)?;

            Ok::<_, BiometalError>((n1, n2, n3))
        })
        .await
        .unwrap()
        .unwrap();

        assert_eq!(n1, 5);
        assert_eq!(n2, 0);
        assert_eq!(n3, 0);
    }

    #[tokio::test]
    async fn test_retry_on_transient_failure() {
        let mock_server = MockServer::start().await;
        let test_data = b"Retry test";

        // Mock HEAD request
        Mock::given(method("HEAD"))
            .and(path("/retry.txt"))
            .respond_with(
                ResponseTemplate::new(200)
                    .insert_header("Content-Length", test_data.len().to_string()),
            )
            .mount(&mock_server)
            .await;

        // First request fails (500)
        Mock::given(method("GET"))
            .and(path("/retry.txt"))
            .respond_with(ResponseTemplate::new(500))
            .up_to_n_times(1)
            .mount(&mock_server)
            .await;

        // Second request succeeds
        Mock::given(method("GET"))
            .and(path("/retry.txt"))
            .respond_with(
                ResponseTemplate::new(206)
                    .set_body_bytes(test_data)
                    .insert_header(
                        "Content-Range",
                        format!("bytes 0-{}/{}", test_data.len() - 1, test_data.len()),
                    ),
            )
            .mount(&mock_server)
            .await;

        let url = format!("{}/retry.txt", mock_server.uri());

        let result = tokio::task::spawn_blocking(move || {
            let client = HttpClient::new()?;
            client.fetch_range(&url, 0, test_data.len() as u64)
        })
        .await
        .unwrap();

        let data = result.unwrap();
        assert_eq!(data.as_ref(), test_data);
    }
}
