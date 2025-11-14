"""
AlphaFold3 Server Client (No Installation Required!)

Uses the free AlphaFold3 Server API - just need an API key.
Free tier: 20 predictions per day (perfect for pilot studies)

Get your API key at: https://golgi.sandbox.google.com/
"""

import logging
import json
import time
import requests
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class AF3ServerConfig:
    """Configuration for AlphaFold3 Server"""
    api_key: str
    base_url: str = "https://alphafoldserver.com/api/v1"
    max_retries: int = 3
    poll_interval: int = 60  # seconds


class AlphaFold3ServerClient:
    """
    Easy-to-use client for AlphaFold3 Server

    No installation required - just an API key!
    """

    def __init__(self, api_key: str):
        """
        Initialize client with API key

        Args:
            api_key: Your AlphaFold3 Server API key
                    Get one at: https://golgi.sandbox.google.com/
        """
        self.api_key = api_key
        self.base_url = "https://alphafoldserver.com/api/v1"
        self.session = requests.Session()
        self.session.headers.update({
            'Authorization': f'Bearer {api_key}',
            'Content-Type': 'application/json'
        })

    def submit_prediction(self, input_json: Path) -> str:
        """
        Submit a prediction job

        Args:
            input_json: Path to AlphaFold3 JSON input file

        Returns:
            Job ID
        """
        logger.info(f"Submitting prediction: {input_json.name}")

        # Load input
        with open(input_json) as f:
            input_data = json.load(f)

        # Submit to API
        url = f"{self.base_url}/fold"

        try:
            response = self.session.post(url, json=input_data, timeout=30)
            response.raise_for_status()

            result = response.json()
            job_id = result.get('id')

            logger.info(f"✓ Job submitted: {job_id}")
            return job_id

        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to submit job: {e}")
            raise

    def get_job_status(self, job_id: str) -> Dict:
        """
        Check job status

        Returns:
            Dictionary with status, progress, etc.
        """
        url = f"{self.base_url}/fold/{job_id}"

        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to get status: {e}")
            raise

    def wait_for_completion(
        self,
        job_id: str,
        max_wait: int = 3600,
        poll_interval: int = 60
    ) -> Dict:
        """
        Wait for job to complete

        Args:
            job_id: Job identifier
            max_wait: Maximum wait time in seconds
            poll_interval: How often to check (seconds)

        Returns:
            Final job result
        """
        logger.info(f"Waiting for job {job_id}...")
        start_time = time.time()

        while time.time() - start_time < max_wait:
            status = self.get_job_status(job_id)

            state = status.get('state')

            if state == 'COMPLETED':
                logger.info(f"✓ Job {job_id} completed!")
                return status

            elif state == 'FAILED':
                error = status.get('error', 'Unknown error')
                logger.error(f"✗ Job {job_id} failed: {error}")
                raise RuntimeError(f"Job failed: {error}")

            elif state in ['PENDING', 'RUNNING']:
                elapsed = int(time.time() - start_time)
                logger.info(f"Job {job_id} {state.lower()} - elapsed: {elapsed}s")
                time.sleep(poll_interval)

            else:
                logger.warning(f"Unknown state: {state}")
                time.sleep(poll_interval)

        raise TimeoutError(f"Job {job_id} did not complete within {max_wait}s")

    def download_results(self, job_id: str, output_dir: Path):
        """
        Download prediction results

        Args:
            job_id: Job identifier
            output_dir: Where to save results
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get job details
        status = self.get_job_status(job_id)

        if status.get('state') != 'COMPLETED':
            raise RuntimeError(f"Job not completed: {status.get('state')}")

        # Download files
        files = status.get('files', [])

        for file_info in files:
            file_url = file_info.get('url')
            file_name = file_info.get('name')

            if not file_url or not file_name:
                continue

            logger.info(f"Downloading {file_name}...")

            try:
                response = self.session.get(file_url, timeout=120)
                response.raise_for_status()

                output_file = output_dir / file_name
                with open(output_file, 'wb') as f:
                    f.write(response.content)

                logger.info(f"✓ Saved to {output_file}")

            except Exception as e:
                logger.error(f"Failed to download {file_name}: {e}")

    def predict(
        self,
        input_json: Path,
        output_dir: Path,
        wait: bool = True
    ) -> Dict:
        """
        Complete prediction workflow: submit + wait + download

        Args:
            input_json: AlphaFold3 JSON input
            output_dir: Where to save results
            wait: Whether to wait for completion

        Returns:
            Job result dictionary
        """
        # Submit
        job_id = self.submit_prediction(input_json)

        if not wait:
            logger.info(f"Job submitted: {job_id} (not waiting)")
            return {'job_id': job_id, 'state': 'PENDING'}

        # Wait for completion
        result = self.wait_for_completion(job_id)

        # Download results
        self.download_results(job_id, output_dir)

        return result

    def predict_batch(
        self,
        input_files: List[Path],
        output_dir: Path,
        max_parallel: int = 20,  # Free tier limit
        wait_between: int = 2  # Avoid rate limiting
    ) -> Dict[str, Dict]:
        """
        Batch prediction (respects free tier limits)

        Args:
            input_files: List of JSON input files
            output_dir: Base output directory
            max_parallel: Maximum concurrent jobs (20 for free tier)
            wait_between: Seconds to wait between submissions

        Returns:
            Dictionary mapping filename to results
        """
        results = {}
        pending_jobs = {}

        logger.info(f"Submitting {len(input_files)} predictions...")
        logger.info(f"Free tier limit: {max_parallel} jobs/day")

        # Submit all jobs (up to daily limit)
        for i, input_file in enumerate(input_files[:max_parallel], 1):
            try:
                logger.info(f"Submitting {i}/{min(len(input_files), max_parallel)}: {input_file.name}")
                job_id = self.submit_prediction(input_file)
                pending_jobs[job_id] = {
                    'input_file': input_file,
                    'complex_id': input_file.stem
                }

                # Rate limiting
                if i < len(input_files):
                    time.sleep(wait_between)

            except Exception as e:
                logger.error(f"Failed to submit {input_file.name}: {e}")
                results[input_file.stem] = {'error': str(e)}

        logger.info(f"\n✓ Submitted {len(pending_jobs)} jobs")
        logger.info("Waiting for completion...")

        # Wait for all jobs to complete
        completed = 0
        while pending_jobs:
            for job_id in list(pending_jobs.keys()):
                try:
                    status = self.get_job_status(job_id)
                    state = status.get('state')

                    if state == 'COMPLETED':
                        job_info = pending_jobs[job_id]
                        complex_id = job_info['complex_id']

                        # Download results
                        job_output_dir = output_dir / complex_id
                        self.download_results(job_id, job_output_dir)

                        results[complex_id] = {
                            'job_id': job_id,
                            'status': 'completed',
                            'output_dir': str(job_output_dir)
                        }

                        del pending_jobs[job_id]
                        completed += 1

                        logger.info(f"✓ Completed {completed}/{len(results)} predictions")

                    elif state == 'FAILED':
                        job_info = pending_jobs[job_id]
                        complex_id = job_info['complex_id']
                        error = status.get('error', 'Unknown error')

                        results[complex_id] = {
                            'job_id': job_id,
                            'status': 'failed',
                            'error': error
                        }

                        del pending_jobs[job_id]
                        logger.error(f"✗ Job {job_id} failed: {error}")

                except Exception as e:
                    logger.warning(f"Error checking job {job_id}: {e}")

            if pending_jobs:
                time.sleep(60)  # Check every minute

        logger.info(f"\n✓ All jobs complete: {completed} successful")
        return results


def main():
    """
    Example usage and API key setup
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="AlphaFold3 Server Client (No Installation Required!)"
    )
    parser.add_argument(
        '--api-key',
        type=str,
        help='Your AlphaFold3 Server API key'
    )
    parser.add_argument(
        '--input',
        type=Path,
        help='Input JSON file or directory'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('results/predictions'),
        help='Output directory'
    )
    parser.add_argument(
        '--setup',
        action='store_true',
        help='Show setup instructions'
    )

    args = parser.parse_args()

    if args.setup or not args.api_key:
        print("""
╔════════════════════════════════════════════════════════════════╗
║         AlphaFold3 Server Setup (Zero Installation!)          ║
╚════════════════════════════════════════════════════════════════╝

Step 1: Get your FREE API key
   → Go to: https://golgi.sandbox.google.com/
   → Sign in with Google account
   → Create API key
   → Copy the key

Step 2: Save your API key
   → Create file: config/af3_api_key.txt
   → Paste your API key in the file

   Or set environment variable:
   export AF3_API_KEY="your-api-key-here"

Step 3: Run predictions!
   python -m src.prediction.alphafold_server_client \\
       --api-key $(cat config/af3_api_key.txt) \\
       --input data/processed/structures/af3_inputs \\
       --output-dir results/predictions

Free Tier Limits:
   • 20 predictions per day
   • Perfect for pilot studies
   • No GPU or installation needed!
   • Results in ~5-15 minutes per prediction

Need more predictions?
   • Paid tier: Unlimited predictions
   • Or run locally with GPU (see README.md)

Questions? Check the README.md or API docs.
        """)
        return

    # Initialize client
    client = AlphaFold3ServerClient(api_key=args.api_key)

    # Check if input is file or directory
    input_path = Path(args.input)

    if input_path.is_file():
        # Single prediction
        result = client.predict(input_path, args.output_dir)
        print(f"\n✓ Prediction complete: {result.get('state')}")

    elif input_path.is_dir():
        # Batch predictions
        input_files = list(input_path.glob("*.json"))
        print(f"\nFound {len(input_files)} input files")
        print("Free tier limit: 20 per day")

        confirm = input(f"\nSubmit {min(len(input_files), 20)} predictions? [y/N]: ")
        if confirm.lower() != 'y':
            print("Cancelled")
            return

        results = client.predict_batch(input_files, args.output_dir)

        # Summary
        successful = sum(1 for r in results.values() if r.get('status') == 'completed')
        print(f"\n✓ Batch complete: {successful}/{len(results)} successful")

    else:
        print(f"Error: {input_path} not found")


if __name__ == "__main__":
    main()
